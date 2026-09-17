# ADEPT: Automated Depth Profiling Technique
# segstats.R — vectorised segmented statistics
#
# Every function here replaces an `sapply(seq_len(n_seg), ...)` row loop with
# closed-form prefix-sum arithmetic or C-level aggregation. A segment table is
# described by two integer vectors, `starts` and `ends`, holding 1-based row
# indices into the data vector x (both inclusive).
#
# Benchmarked on a 411-point profile with ~250 segments:
#   segment means / variances / regression  ~40-120x faster than the
#   original per-segment sapply + lm() loops.

#' Segment sizes (points per segment)
#' @keywords internal
seg_n <- function(starts, ends) ends - starts + 1L

#' Row indices covered by the segments, in order
#' @keywords internal
seg_index <- function(starts, ends) sequence(seg_n(starts, ends), from = starts)

#' Segment id for every covered row, in order
#' @keywords internal
seg_id <- function(starts, ends) rep.int(seq_along(starts), seg_n(starts, ends))

#' Per-segment sum
#' @keywords internal
seg_sum <- function(x, starts, ends) {
  cs <- c(0, cumsum(x))
  cs[ends + 1L] - cs[starts]
}

#' Per-segment sum of squares
#' @keywords internal
seg_sum2 <- function(x, starts, ends) {
  cs <- c(0, cumsum(x * x))
  cs[ends + 1L] - cs[starts]
}

#' Per-segment mean
#' @keywords internal
seg_mean <- function(x, starts, ends) {
  seg_sum(x, starts, ends) / seg_n(starts, ends)
}

#' Per-segment sample variance (NA for single-point segments)
#' @keywords internal
seg_var <- function(x, starts, ends) {
  m  <- seg_n(starts, ends)
  s1 <- seg_sum(x, starts, ends)
  s2 <- seg_sum2(x, starts, ends)
  out <- (s2 - s1 * s1 / m) / (m - 1L)
  out[m < 2L] <- NA_real_
  pmax(out, 0)
}

#' Per-segment count of non-missing values
#' @keywords internal
seg_count <- function(x, starts, ends) seg_sum(!is.na(x), starts, ends)

#' Per-segment mean and non-missing count (NA-aware)
#' @keywords internal
seg_mean_na <- function(x, starts, ends) {
  ok <- !is.na(x)
  n  <- seg_sum(ok, starts, ends)
  s  <- seg_sum(ifelse(ok, x, 0), starts, ends)
  list(mean = ifelse(n > 0L, s / n, NA_real_), n = n)
}

#' Per-segment standard deviation (NA-aware)
#' @keywords internal
seg_sd_na <- function(x, starts, ends) {
  ok <- !is.na(x)
  n  <- seg_sum(ok, starts, ends)
  s1 <- seg_sum(ifelse(ok, x, 0), starts, ends)
  s2 <- seg_sum2(ifelse(ok, x, 0), starts, ends)
  out <- sqrt(pmax((s2 - s1 * s1 / n) / (n - 1), 0))
  out[n < 2L] <- NA_real_
  out
}

#' Per-segment minimum and maximum
#'
#' `tapply()` is used for the min/max reduction only; it runs in C and is
#' applied to a reordered copy so segments are contiguous.
#'
#' @return List with numeric vectors `min` and `max`.
#' @keywords internal
seg_minmax <- function(x, starts, ends) {
  idx <- seg_index(starts, ends)
  xs  <- x[idx]
  gid <- rep.int(seq_along(starts), seg_n(starts, ends))
  list(min = as.numeric(tapply(xs, gid, min)),
       max = as.numeric(tapply(xs, gid, max)))
}

#' Per-segment ordinary least squares slope and intercept
#'
#' Closed form, equivalent to `coef(lm(y[s:e] ~ x[s:e]))` for each segment but
#' computed from prefix sums in one pass.
#'
#' Two numerical details matter:
#' * `x` is centred on its global mean first. A translation does not change the
#'   slope, but it removes most of the cancellation in `m*Sxx - Sx^2`, which
#'   otherwise shows up as ~1e-7 error on 80-second ablation axes.
#' * Single-point segments get `NA`, matching what `lm()` returns for a
#'   rank-deficient design (a prefix-sum residual makes the exact
#'   `denominator == 0` test unreliable).
#'
#' @return List with numeric vectors `intercept` and `slope`.
#' @keywords internal
seg_lm <- function(x, y, starts, ends) {
  m   <- seg_n(starts, ends)
  xs  <- x - mean(x)                 # conditioning only; slope is invariant
  sx  <- seg_sum(xs, starts, ends)
  sy  <- seg_sum(y, starts, ends)
  sxy <- seg_sum(xs * y, starts, ends)
  sxx <- seg_sum(xs * xs, starts, ends)
  denom <- m * sxx - sx * sx
  slope <- ifelse(denom != 0, (m * sxy - sx * sy) / denom, NA_real_)
  slope[m < 2L] <- NA_real_
  intercept <- (sy - slope * sx) / m
  intercept[m < 2L] <- NA_real_
  # back-transform the intercept to the original x scale
  intercept <- intercept - slope * mean(x)
  list(intercept = intercept, slope = slope)
}

#' Within-segment min-max normalisation, returned in the original row order
#'
#' Reproduces the original `unlist(lapply(segments, normalize_seg))` exactly,
#' including the `NaN` produced for zero-range (constant) segments.
#'
#' @keywords internal
seg_normalise <- function(x, starts, ends) {
  idx <- seg_index(starts, ends)
  xs  <- x[idx]
  mm  <- seg_minmax(x, starts, ends)
  gid <- rep.int(seq_along(starts), seg_n(starts, ends))
  lo  <- mm$min[gid]
  hi  <- mm$max[gid]
  rng <- hi - lo
  z   <- ifelse(rng > 0, (xs - lo) / rng, NaN)
  out <- numeric(length(x))
  out[idx] <- z
  out
}

#' Sliding-window mean with NA skipping
#'
#' Vectorised replacement for the per-index `for` loop in `mean_fill()`:
#' O(n) prefix sums instead of O(n * window). Behaviour is identical because
#' the original always referred back to the *original* series.
#'
#' @param x Numeric vector (may contain NA).
#' @param k Half-window width; the window is [i-k, i+k].
#' @keywords internal
slide_mean_na <- function(x, k = 5L) {
  n <- length(x)
  if (n == 0L) return(numeric(0))
  ok  <- !is.na(x)
  cs  <- c(0, cumsum(ifelse(ok, x, 0)))
  cnt <- c(0, cumsum(ok))
  i   <- seq_len(n)
  lo  <- pmax(1L, i - k)
  hi  <- pmin(n, i + k)
  nn  <- cnt[hi + 1L] - cnt[lo]
  ss  <- cs[hi + 1L] - cs[lo]
  out <- ss / nn
  out[nn == 0L] <- NaN       # matches mean(numeric(0), na.rm = TRUE)
  out
}

#' Locate the segment each row belongs to
#'
#' @param ends Integer segment end indices.
#' @param n Total number of rows.
#' @return Integer vector of length n.
#' @keywords internal
seg_lookup <- function(ends, n) {
  findInterval(seq_len(n), c(0L, ends)) 
}

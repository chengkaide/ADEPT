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

#' Per-segment inverse-variance weighted mean, standard error and MSWD
#'
#' Used when the input carries a per-point 1-sigma. The plateau age is then the
#' inverse-variance weighted mean of the measured ages rather than the mean of
#' the smoothed curve, because the weights already hold the information the
#' LOESS fit was there to supply.
#'
#' MSWD (mean square of weighted deviates) asks whether the scatter of the
#' points is consistent with the errors they were given:
#'
#'     MSWD = sum(w_i * (x_i - xbar)^2) / (n - 1),   w_i = 1 / sigma_i^2
#'
#' Around 1 it is: the spread matches the errors and the plateau can be read as
#' a single age. Much above 1 means either real age heterogeneity inside the
#' plateau or sigmas that are too small (the usual cause is counting statistics
#' quoted without the external reproducibility). Much below 1 means sigmas that
#' are too large.
#'
#' `prob` is the probability of seeing a value at least this large under that
#' consistency assumption; below about 0.05 is the conventional flag.
#'
#' Points with a non-finite age or a non-positive sigma are dropped rather than
#' allowed to fail the whole segment, matching the counterpart implementation
#' in the Python reduction.
#'
#' @return List of numeric vectors `mean`, `se`, `mswd`, `prob` and the integer
#'   `n` of usable points per segment.
#' @keywords internal
seg_weighted_mean <- function(x, s, starts, ends) {
  ok <- is.finite(x) & is.finite(s) & s > 0
  w   <- ifelse(ok, 1 / s^2, 0)
  wx  <- ifelse(ok, x / s^2, 0)
  wxx <- ifelse(ok, x * x / s^2, 0)

  W   <- seg_sum(w,   starts, ends)
  WX  <- seg_sum(wx,  starts, ends)
  WXX <- seg_sum(wxx, starts, ends)
  n   <- seg_sum(ok,  starts, ends)

  mu <- WX / W
  se <- 1 / sqrt(W)
  # Weighted residual sum of squares, in the one-pass form. Safe here because W
  # is a sum of positive weights and cannot be near zero for a usable segment;
  # points with sigma <= 0 were already given zero weight.
  ss   <- WXX - WX * WX / W
  mswd <- ss / (n - 1)
  prob <- stats::pf(mswd, n - 1, Inf, lower.tail = FALSE)

  # A single usable point has no scatter to test, so MSWD and the standard
  # error are undefined. The weighted mean itself is still well defined and is
  # kept, matching how a one-point segment behaves when no sigmas are supplied.
  bad <- n < 2L
  mu[n == 0L] <- NA_real_
  se[bad]     <- NA_real_
  mswd[bad]   <- NA_real_
  prob[bad]   <- NA_real_

  list(mean = mu, se = se, mswd = mswd, prob = prob, n = as.integer(n))
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

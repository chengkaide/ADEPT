# ADEPT: Automated Depth Profiling Technique
# filtering.R — Four-step plateau filtering
#
# Internal functions implementing the multi-step filtering cascade.

#' Step 1: Filter by age range
#'
#' Keeps plateaus whose Segment_Mean falls within [min_age, max_age].
#'
#' @param segments Segment data.frame
#' @param min_age Minimum valid age (Ma)
#' @param max_age Maximum valid age (Ma)
#' @return Numeric vector (Filter_1): age value or NA
#' @keywords internal
filter_age_range <- function(segments, min_age = 0, max_age = 4540) {
  ifelse(as.numeric(segments$Segment_Mean) <= max_age &
         as.numeric(segments$Segment_Mean) >= min_age,
         segments$Segment_Mean, NA)
}

#' Step 2: Filter by plateau variance
#'
#' @param segments Segment data.frame
#' @param filter_1 Result from filter_age_range()
#' @param threshold Variance threshold (default 0.1192)
#' @return Numeric vector (Filter_2): age value or NA
#' @keywords internal
filter_variance <- function(segments, filter_1, threshold = 0.1192) {
  ifelse(as.numeric(segments$Variance) <= threshold, filter_1, NA)
}

#' Step 3: Filter by minimum plateau resolution
#'
#' Removes plateaus shorter than the minimum resolution (default >= 5s).
#'
#' @param segments Segment data.frame
#' @param filter_2 Result from filter_variance()
#' @param min_res Minimum plateau duration in seconds (NULL or < 5 defaults to 5)
#' @return Numeric vector (Filter_3): age value or NA
#' @keywords internal
filter_resolution <- function(segments, filter_2, min_res = NULL,
                              skip_first = TRUE) {
  if (is.null(min_res) || length(min_res) != 1 || is.na(min_res) ||
      min_res < 5) min_res <- 5

  result <- filter_2
  n <- nrow(segments)
  if (n <= 1) return(result)   # nothing to compare / single segment

  # NOTE: the first segment is exempt from the minimum-duration rule.
  # It carries the onset of the ablation signal, which is often truncated
  # by `lower_ablation_time` and would otherwise always be discarded.
  idx <- if (isTRUE(skip_first)) 2:n else seq_len(n)
  result[idx] <- ifelse(as.numeric(segments$Time_step[idx]) >= min_res,
                        filter_2[idx], NA)
  result
}

#' Longest monotonic subsequence (kept indices)
#'
#' Dynamic programming over the plateau ages. Used to keep the largest set of
#' plateaus that is consistent with the expected age trend, instead of
#' collapsing to two hand-picked plateaus.
#'
#' @param a Numeric vector of plateau ages.
#' @param increasing TRUE for a forward (age increasing) trend.
#' @param tol Allowed relative reversal; a step of up to \code{tol * a[i]}
#'   against the trend still counts as monotonic.
#' @return Integer vector of indices, in ascending order.
#' @keywords internal
longest_monotonic <- function(a, increasing = TRUE, tol = 0) {
  n <- length(a)
  if (n == 0L) return(integer(0))
  if (n == 1L) return(1L)

  ok <- function(i, j) {
    # i comes before j
    if (increasing) a[j] >= a[i] * (1 - tol) else a[j] <= a[i] * (1 + tol)
  }

  dp   <- rep(1L, n)
  prev <- rep(0L, n)
  for (j in seq_len(n)[-1]) {
    for (i in seq_len(j - 1L)) {
      if (ok(i, j) && dp[i] + 1L > dp[j]) {
        dp[j]   <- dp[i] + 1L
        prev[j] <- i
      }
    }
  }
  k <- which.max(dp)
  out <- integer(0)
  while (k > 0L) {
    out <- c(k, out)
    k <- prev[k]
  }
  out
}

#' Step 4: Directional plateau selection (Forward / Reverse)
#'
#' Forward: keep the largest set of plateaus whose ages do not decrease.
#' Reverse: keep the largest set whose ages do not increase.
#'
#' Two methods are available:
#'   "monotonic" (default) — longest monotonic subsequence, with a tolerance
#'     so that sub-analytical reversals do not break the trend. This keeps
#'     every plateau that is consistent with the trend.
#'   "strict" — the original v1.1.0 behaviour: require a perfectly monotonic
#'     sequence, otherwise fall back to the first plateau plus the one with
#'     the smallest variance (at most two plateaus). Kept for reproducing
#'     earlier results.
#'
#' @param segments Segment data.frame
#' @param filter_3 Result from filter_resolution()
#' @param direction "Forward" (keep ascending) or "Reverse" (keep descending)
#' @param method "monotonic" (default) or "strict"
#' @param tol Relative reversal tolerated by the "monotonic" method
#' @return Numeric vector (Filter_4): final age value or NA
#' @keywords internal
filter_direction <- function(segments, filter_3,
                             direction = c("Forward", "Reverse"),
                             method = c("monotonic", "strict"),
                             tol = 0.02) {
  direction <- match.arg(direction)
  method    <- match.arg(method)
  filter_4  <- rep(NA_real_, nrow(segments))

  keep_idx <- which(!is.na(filter_3))
  if (length(keep_idx) == 0L) return(filter_4)
  vals <- as.numeric(filter_3[keep_idx])
  increasing <- identical(direction, "Forward")

  if (identical(method, "strict")) {
    # ---- original v1.1.0 behaviour -----------------------------------------
    is_ascending  <- all(diff(vals) >= 0)
    is_descending <- all(diff(vals) < 0)
    if ((increasing && is_ascending) || (!increasing && is_descending)) {
      filter_4[keep_idx] <- vals
      return(filter_4)
    }
    filter_4[keep_idx[1]] <- vals[1]
    other <- if (increasing) !is_descending else !is_ascending
    if (other) {
      remaining <- keep_idx[-1]
      if (length(remaining) > 0L) {
        best <- remaining[which.min(segments$Variance[remaining])]
        filter_4[best] <- filter_3[best]
      }
    }
    return(filter_4)
  }

  # ---- v1.2.0: longest monotonic subsequence -------------------------------
  sel <- longest_monotonic(vals, increasing = increasing, tol = tol)
  filter_4[keep_idx[sel]] <- vals[sel]
  filter_4
}

#' Apply all four filtering steps
#'
#' @param segments Segment data.frame
#' @param min_age, max_age Age limits (Ma)
#' @param var_threshold Variance threshold
#' @param min_res Minimum plateau resolution (seconds)
#' @param direction "Forward" or "Reverse"
#' @return Updated segments data.frame with Filter_1..Filter_4 and derived columns
#' @keywords internal
apply_filters <- function(segments, min_age = 0, max_age = 4540,
                          var_threshold = 0.1192, min_res = NULL,
                          direction = c("Forward", "Reverse"),
                          direction_method = c("monotonic", "strict"),
                          direction_tolerance = 0.02) {
  direction        <- match.arg(direction)
  direction_method <- match.arg(direction_method)

  segments$Filter_1 <- filter_age_range(segments, min_age, max_age)
  segments$Filter_2 <- filter_variance(segments, segments$Filter_1, var_threshold)
  segments$Filter_3 <- filter_resolution(segments, segments$Filter_2, min_res)
  segments$Filter_4 <- filter_direction(segments, segments$Filter_3, direction,
                                        method = direction_method,
                                        tol = direction_tolerance)

  segments$Final_total_uncertainty <- ifelse(!is.na(segments$Filter_4),
                                             segments$Total_uncertainty, NA)
  # v1.2.0: the uncertainty that includes the decay constants
  segments$Final_total_uncertainty_full <-
    ifelse(!is.na(segments$Filter_4), segments$Total_uncertainty_full, NA)

  # v1.2.0 rename: the old name was misleading, it is the total number of
  # segments produced by PELT, not the number of confirmed plateaus.
  segments$Total_Segments  <- nrow(segments)
  segments$Plateau_Numbers <- nrow(segments)   # kept for backwards compat

  segments$Final_Serial_Number <- NA
  non_na <- which(!is.na(segments$Filter_4))
  segments$Final_Serial_Number[non_na] <- seq_along(non_na)
  segments$Final_Age <- segments$Filter_4
  segments$Final_step_Numbers <- sum(!is.na(segments$Filter_4))
  segments$Confirmed_Plateaus <- sum(!is.na(segments$Filter_4))

  return(segments)
}

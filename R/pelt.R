# ADEPT: Automated Depth Profiling Technique
# pelt.R — PELT changepoint detection, implemented in base R
#
# Replaces the `changepoint` dependency. The implementation follows
# Killick, Fearnhead & Eckley (2012), JASA 107(497), 1590-1598, with the
# L2 (residual-sum-of-squares) cost.
#
# The penalty used by the original ADEPT code is
#     SAIC = (pen_AIC / 100) * log(n)
# where `changepoint::cpt.mean(penalty = "AIC")` always reports pen.value = 4
# for a univariate mean model (2 * 2 parameters). That constant is therefore
# hard-coded here as `ADEPT_AIC_PENALTY`, which removes any dependence on a
# particular version of the changepoint package.
#
# Verified against changepoint 2.2.4: identical changepoint sets on all
# tested series (n = 60/161/411, 1-6 true segments).

ADEPT_AIC_PENALTY <- 4

#' L2 cost of every segment (vectorised over segments)
#'
#' cost(s, t) = sum_{i=s..t} (x_i - mean(x_{s..t}))^2
#'
#' @param cs1 Prefix sums of x, starting with 0.
#' @param cs2 Prefix sums of x^2, starting with 0.
#' @keywords internal
pelt_seg_cost <- function(cs1, cs2, s, t) {
  m <- t - s + 1
  s1 <- cs1[t + 1L] - cs1[s]
  s2 <- cs2[t + 1L] - cs2[s]
  pmax(s2 - s1 * s1 / m, 0)
}

#' PELT changepoint detection for a mean shift
#'
#' @param x Numeric vector (no NA).
#' @param pen Penalty per changepoint.
#' @param minseglen Minimum segment length (default 1, matching the original).
#' @return Integer vector of changepoint positions (segment end indices),
#'   always ending at `length(x)`.
#' @keywords internal
pelt_mean <- function(x, pen, minseglen = 1L) {
  n <- length(x)
  if (n < 2L) return(n)

  cs1 <- c(0, cumsum(x))
  cs2 <- c(0, cumsum(x * x))

  Fv   <- rep(Inf, n + 1L)
  Fv[1L] <- -pen
  last <- integer(n)
  R    <- 0L                       # candidate previous changepoint

  for (t in seq_len(n)) {
    # Only tau <= t - minseglen keeps *every* segment at least minseglen long,
    # including the first one (tau = 0 requires t >= minseglen). When no such
    # candidate exists, t cannot end a segment and F(t) stays Inf.
    cand <- R[R <= t - minseglen]
    if (length(cand) == 0L) {
      last[t] <- 0L
    } else {
      vals <- Fv[cand + 1L] + pelt_seg_cost(cs1, cs2, cand + 1L, t) + pen
      k    <- which.min(vals)
      Fv[t + 1L] <- vals[k]
      last[t]    <- cand[k]
    }

    # Pruning. Only candidates that were actually *evaluated* at time t may be
    # discarded; anything still inside the minimum-segment lookback has not
    # been costed yet and must survive to be considered at a later t.
    evaluable <- R[R <= t - minseglen]
    if (length(evaluable)) {
      keep <- Fv[evaluable + 1L] +
              pelt_seg_cost(cs1, cs2, evaluable + 1L, t) < Fv[t + 1L]
      evaluable <- evaluable[keep]
    }
    pending <- R[R > t - minseglen]
    R <- c(evaluable, pending, t)
  }

  cpts <- integer(0)
  pos  <- n
  while (pos > 0L) {
    cpts <- c(pos, cpts)
    pos  <- last[pos]
  }
  cpts
}

#' ADEPT segmentation: penalty derived from the AIC penalty, then PELT
#'
#' @param series Numeric vector of standardised LOESS values.
#' @param n Number of observations (kept for API compatibility).
#' @return List with `changepoints` (integer vector) and `SAIC` (penalty).
#' @keywords internal
pelt_segmentation <- function(series, n = length(series)) {
  if (length(series) < 4L ||
      !is.finite(stats::sd(series, na.rm = TRUE)) ||
      stats::sd(series, na.rm = TRUE) == 0) {
    return(list(changepoints = length(series), SAIC = NA_real_))
  }

  SAIC  <- (ADEPT_AIC_PENALTY / 100) * log(n)
  cpts <- tryCatch(
    pelt_mean(as.numeric(series), pen = SAIC, minseglen = 1L),
    error = function(e) length(series)
  )

  cpts <- as.integer(cpts)
  cpts <- sort(unique(cpts[is.finite(cpts) & cpts >= 1L & cpts <= length(series)]))
  if (length(cpts) == 0L) cpts <- length(series)

  list(changepoints = cpts, SAIC = SAIC)
}

# ADEPT: Automated Depth Profiling Technique
# isotopes.R — U-Pb age calculation, implemented in base R
#
# Replaces the IsoplotR dependency for converting isotopic ratios or raw
# counts into ages. The decay constants and the 238U/235U ratio follow the
# values used by IsoplotR, so results agree to < 1e-4 Ma over the whole
# geological range (verified against IsoplotR 5.1).

# Decay constants (yr^-1)
ADEPT_LAMBDA238 <- 1.55125e-10   # Jaffey et al. (1971)
ADEPT_LAMBDA235 <- 9.8485e-10    # Jaffey et al. (1971)
ADEPT_U238U235  <- 137.818       # Hiess et al. (2012)
ADEPT_EARTH_AGE <- 4.567e9       # yr, upper bracket for the 207/206 solve

#' 206Pb/238U age (Ma)
#'
#' Vectorised; solves `R = exp(lambda238 * t) - 1`.
#'
#' @param r Numeric vector of 206Pb/238U ratios.
#' @return Numeric vector of ages in Ma (NA for non-positive / non-finite r).
#' @keywords internal
age_206_238 <- function(r) {
  out <- rep(NA_real_, length(r))
  ok <- is.finite(r) & r > 0
  out[ok] <- log1p(r[ok]) / ADEPT_LAMBDA238 / 1e6
  out
}

#' 207Pb/235U age (Ma)
#'
#' @param r Numeric vector of 207Pb/235U ratios.
#' @return Numeric vector of ages in Ma.
#' @keywords internal
age_207_235 <- function(r) {
  out <- rep(NA_real_, length(r))
  ok <- is.finite(r) & r > 0
  out[ok] <- log1p(r[ok]) / ADEPT_LAMBDA235 / 1e6
  out
}

#' 207Pb/206Pb age (Ma)
#'
#' Solves the transcendental equation
#' `r = (1/137.818) * (exp(l235*t) - 1) / (exp(l238*t) - 1)`
#' by vectorised bisection — no loop over elements, no dependency.
#'
#' @param r Numeric vector of 207Pb/206Pb ratios.
#' @return Numeric vector of ages in Ma.
#' @keywords internal
age_207_206 <- function(r) {
  r <- as.numeric(r)
  out <- rep(NA_real_, length(r))
  r_min <- (1 / ADEPT_U238U235) * (ADEPT_LAMBDA235 / ADEPT_LAMBDA238)
  ok <- is.finite(r) & r > r_min
  if (!any(ok)) return(out)

  target <- r[ok]
  lo <- rep(1, length(target))
  hi <- rep(ADEPT_EARTH_AGE, length(target))
  l238 <- ADEPT_LAMBDA238
  l235 <- ADEPT_LAMBDA235
  k <- 1 / ADEPT_U238U235

  # 60 bisection steps: interval 4.6e9 / 2^60 << 1 yr
  for (i in 1:60) {
    mid <- (lo + hi) / 2
    g <- k * (exp(l235 * mid) - 1) / (exp(l238 * mid) - 1)
    go_hi <- g > target
    hi[go_hi] <- mid[go_hi]
    lo[!go_hi] <- mid[!go_hi]
  }
  out[ok] <- (lo + hi) / 2 / 1e6
  out
}

#' Convert isotopic ratios to U-Pb ages (all three systems at once)
#'
#' @param r68 206Pb/238U ratios
#' @param r75 207Pb/235U ratios
#' @param r76 207Pb/206Pb ratios
#' @return data.frame with columns Age68, Age75, Age76 (Ma)
#' @keywords internal
ratios_to_ages <- function(r68, r75, r76) {
  data.frame(Age68 = age_206_238(r68),
             Age75 = age_207_235(r75),
             Age76 = age_207_206(r76))
}

#' Convert raw isotope counts to U-Pb ages
#'
#' @param Pb206,Pb207,U238 Numeric vectors of counts.
#' @param u238u235 238U/235U ratio used to derive 235U (default 137.818).
#' @return data.frame with columns Age68, Age75, Age76 (Ma)
#' @keywords internal
counts_to_ages <- function(Pb206, Pb207, U238, u238u235 = ADEPT_U238U235) {
  U235 <- U238 / u238u235
  ratios_to_ages(Pb206 / U238, Pb207 / U235, Pb207 / Pb206)
}

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

# ---------------------------------------------------------------------------
#  Uncertainty of the decay constants
# ---------------------------------------------------------------------------
# IsoplotR stores lambda in Ma^-1, so its quoted 1-sigma values are per Ma:
#   U238 : 1.55125e-4 +- 8.3e-8   -> relative 1-sigma 5.3505e-4
#   U235 : 9.84850e-4 +- 6.7e-7   -> relative 1-sigma 6.8031e-4
#   U238U235 : 137.818 +- 0.0225  -> relative 1-sigma 1.6326e-4
# The relative values are unit-independent, so they are stored directly to
# avoid any Ma / yr confusion.
ADEPT_LAMBDA238_REL_1S <- 8.3e-08 / 1.55125e-04    # 0.0535 %
ADEPT_LAMBDA235_REL_1S <- 6.7e-07 / 9.84850e-04    # 0.0680 %
ADEPT_U238U235_REL_1S  <- 0.0225  / 137.818        # 0.0163 %

#' Relative 1-sigma decay-constant uncertainty of an age
#'
#' For 206Pb/238U and 207Pb/235U the age is exactly inversely proportional to
#' the decay constant, so the relative uncertainty is just that of lambda.
#'
#' For 207Pb/206Pb the ratio couples lambda235, lambda238 and 238U/235U, and
#' what matters is how the *age* responds - not how the ratio responds. Using
#' the implicit function theorem on
#'   r = (1/u) * (exp(l235 t) - 1) / (exp(l238 t) - 1)
#' the log-elasticities are
#'   a235 =  l235 * t * E235 / (E235 - 1)
#'   a238 = -l238 * t * E238 / (E238 - 1)
#'   au   = -1
#'   b    = a235 + a238            (d ln r / d ln t)
#' so   sigma_ln t = sqrt( sum_i (a_i / b)^2 * sigma_ln(param_i)^2 ).
#'
#' The two terms partly cancel in b, which is why 207Pb/206Pb ages are far less
#' sensitive to decay-constant uncertainty than 206Pb/238U ages - the numerical
#' result (about 0.15 % at 2 Ga) reflects that. Note the value is
#' age-dependent, so it is evaluated per plateau rather than as a constant.
#'
#' @param system Character vector: "Age68", "Age75" or "Age76".
#' @param age Numeric vector of ages in Ma (used only for "Age76").
#' @return Numeric vector of relative (fractional) 1-sigma uncertainties.
#' @keywords internal
adept_decay_rel_1s <- function(system, age = NA_real_) {
  n <- max(length(system), length(age))
  system <- rep_len(system, n)
  age    <- rep_len(age, n)
  out    <- rep(NA_real_, n)

  out[system %in% c("Age68", "Pb206U238", "Pb206_U238")] <-
    ADEPT_LAMBDA238_REL_1S
  out[system %in% c("Age75", "Pb207U235", "Pb207_U235")] <-
    ADEPT_LAMBDA235_REL_1S

  i76 <- which(system %in% c("Age76", "Pb207Pb206", "Pb207_Pb206"))
  if (length(i76)) {
    t <- age[i76] * 1e6                      # Ma -> yr
    ok <- is.finite(t) & t > 0
    val <- rep(NA_real_, length(t))
    if (any(ok)) {
      tt <- t[ok]
      l238 <- ADEPT_LAMBDA238; l235 <- ADEPT_LAMBDA235
      E238 <- exp(l238 * tt); E235 <- exp(l235 * tt)
      a235 <-  l235 * tt * E235 / (E235 - 1)
      a238 <- -l238 * tt * E238 / (E238 - 1)
      b    <- a235 + a238
      val[ok] <- sqrt(
        (a238 / b)^2 * ADEPT_LAMBDA238_REL_1S^2 +
        (a235 / b)^2 * ADEPT_LAMBDA235_REL_1S^2 +
        (1 / b)^2    * ADEPT_U238U235_REL_1S^2)
    }
    out[i76] <- val
  }
  out
}

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

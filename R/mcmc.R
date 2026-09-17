# ADEPT: Automated Depth Profiling Technique
# mcmc.R — Bayesian MCMC posterior analysis (optional)
#
# Uses the `mcp` package to fit a changepoint model and extract
# posterior estimates for each plateau segment.

# Column contract: run_mcmc() ALWAYS returns these columns, even when the
# `mcp` package is missing or the sampler fails. Filling them with NA keeps
# the output schema stable — otherwise `adept()` builds short rows and the
# final `colnames<-` call errors out with
#   "length of 'dimnames' [2] not equal to array extent".
ADEPT_MCMC_COLS <- c("mean", "lower", "upper", "sigma", "Rhat", "n.eff",
                     "mcp_step", "Filter_mcp")

#' Run Bayesian MCMC changepoint analysis
#'
#' Fits an intercept-only changepoint model using the mcp package
#' and returns posterior summaries merged with segment data.
#'
#' @param df data.frame with Time and standardized_loess columns
#' @param segments Segment definitions data.frame
#' @return Updated segments data.frame with MCMC mean, lower, upper, sigma,
#'   Rhat, n.eff and step columns. If `mcp` is unavailable or fitting fails,
#'   the columns are present but NA.
#' @keywords internal
run_mcmc <- function(df, segments) {
  n_segments <- nrow(segments)

  # --- Pre-declare the columns so the schema never changes ---
  for (cn in ADEPT_MCMC_COLS) {
    segments[[cn]] <- rep(NA_real_, n_segments)
  }

  if (!requireNamespace("mcp", quietly = TRUE)) {
    warning("Package 'mcp' is not installed. Skipping MCMC analysis ",
            "(MCMC columns will be NA).\n",
            "  install.packages('mcp')  # requires JAGS: ",
            "https://mcmc-jags.sourceforge.io/", call. = FALSE)
    return(segments)
  }

  if (n_segments < 1) return(segments)

  model <- c(list(stats::as.formula("standardized_loess ~ 1")),
             rep(list(stats::as.formula("~ 1")), n_segments - 1))

  fit <- try(
    mcp::mcp(model, df, par_x = "Time"),
    silent = TRUE
  )
  if (inherits(fit, "try-error")) {
    warning("mcp model construction failed; MCMC columns will be NA. ",
            "Reason: ", conditionMessage(attr(fit, "condition")),
            call. = FALSE)
    return(segments)
  }

  sum_mcp <- try(summary(fit), silent = TRUE)
  if (inherits(sum_mcp, "try-error") || is.null(sum_mcp) || nrow(sum_mcp) == 0) {
    warning("MCMC sampling failed; MCMC columns will be NA.", call. = FALSE)
    return(segments)
  }

  nm <- as.character(sum_mcp$name)
  intercepts_info <- sum_mcp[grep("^int", nm), , drop = FALSE]
  cp_info         <- sum_mcp[grep("^cp",  nm), , drop = FALSE]
  sigma_info      <- sum_mcp[grep("^sigma", nm), , drop = FALSE]

  if (nrow(intercepts_info) == n_segments) {
    segments$mean  <- as.numeric(intercepts_info$mean)
    segments$lower <- as.numeric(intercepts_info$lower)
    segments$upper <- as.numeric(intercepts_info$upper)
    if ("Rhat"   %in% colnames(intercepts_info)) segments$Rhat   <- as.numeric(intercepts_info$Rhat)
    if ("n.eff"  %in% colnames(intercepts_info)) segments$n.eff  <- as.numeric(intercepts_info$n.eff)
    # sigma is a single global parameter; broadcast it to every segment.
    if (nrow(sigma_info) > 0) {
      segments$sigma <- rep(as.numeric(sigma_info$upper[1]), n_segments)
    }
  } else {
    warning("mcp returned ", nrow(intercepts_info), " intercept(s) for ",
            n_segments, " plateau(s); MCMC mean/lower/upper left as NA.",
            call. = FALSE)
  }

  # --- Plateau durations implied by the posterior changepoint locations ---
  mean_vals <- as.numeric(cp_info$mean)
  mean_vals <- mean_vals[is.finite(mean_vals)]
  n_cp <- length(mean_vals)

  if (n_cp >= 1 && n_cp == n_segments - 1) {
    mcp_step <- rep(NA_real_, n_cp + 1)
    t_first <- df$Time[1]
    t_last  <- df$Time[nrow(df)]
    mcp_step[1] <- mean_vals[1] - t_first
    if (n_cp >= 2) mcp_step[2:n_cp] <- diff(mean_vals)
    mcp_step[n_cp + 1] <- t_last - mean_vals[n_cp]
    segments$mcp_step <- mcp_step
  } else if (n_segments == 1) {
    # Single plateau: the whole window is one step.
    segments$mcp_step <- df$Time[nrow(df)] - df$Time[1]
  }

  segments$Filter_mcp <- ifelse(!is.na(segments$Filter_4), segments$mean, NA)

  return(segments)
}

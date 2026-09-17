# ADEPT: Automated Depth Profiling Technique
# sensitivity.R — parameter sensitivity sweep
#
# The plateau count and the resulting ages are sensitive to
# `variance_threshold` and `min_plateau_resolution`. Reviewers ask for that
# sensitivity to be quantified, so it is provided as a first-class function
# rather than as an ad-hoc script.

#' Sweep ADEPT parameters and tabulate the outcome
#'
#' Runs \code{\link{adept}} once per parameter combination and returns a tidy
#' table of the plateau count and ages. Nothing is written to disk.
#'
#' @param file_path Path to the input workbook.
#' @param variance_threshold Numeric vector of variance thresholds to try.
#' @param min_plateau_resolution Numeric vector of minimum plateau durations
#'   (seconds) to try.
#' @param filter_direction Character vector of directions to try.
#' @param calibration_uncertainty Relative 1-sigma reference-material
#'   reproducibility, passed through to \code{\link{adept}}.
#' @param direction_method,direction_tolerance Passed through to
#'   \code{\link{adept}}.
#' @param n_cores Number of workers. 1 (default) keeps everything in the
#'   current process; > 1 uses \pkg{parallel}, which ships with R.
#' @param progress Optional callback \code{function(fraction, detail)}.
#' @param verbose Logical. Print progress?
#' @param ... Further arguments passed to \code{\link{adept}}.
#'
#' @return A data.frame with one row per combination and columns
#'   \code{variance_threshold}, \code{min_plateau_resolution},
#'   \code{filter_direction}, \code{n_plateaus}, \code{ages} (semicolon
#'   separated) and \code{mean_age}.
#'
#' @seealso \code{\link{adept}}
#' @export
#' @examples
#' \dontrun{
#' s <- adept_sensitivity("Input.xlsx")
#' print(s)
#' }
adept_sensitivity <- function(
    file_path,
    variance_threshold     = c(0.05, 0.08, 0.1192, 0.15, 0.25),
    min_plateau_resolution = c(5, 8, 10, 15),
    filter_direction       = c("Forward", "Reverse"),
    calibration_uncertainty = 0.03,
    direction_method       = c("monotonic", "strict"),
    direction_tolerance    = 0.02,
    n_cores                = 1L,
    progress               = NULL,
    verbose                = TRUE,
    ...) {
  direction_method <- match.arg(direction_method)
  grid <- expand.grid(
    variance_threshold     = variance_threshold,
    min_plateau_resolution = min_plateau_resolution,
    filter_direction       = filter_direction,
    stringsAsFactors       = FALSE
  )
  n <- nrow(grid)

  ablate <- function(k) {
    g <- grid[k, ]
    res <- try(suppressWarnings(
      adept(file_path,
            variance_threshold      = g$variance_threshold,
            min_plateau_resolution  = g$min_plateau_resolution,
            filter_direction        = g$filter_direction,
            calibration_uncertainty = calibration_uncertainty,
            direction_method        = direction_method,
            direction_tolerance     = direction_tolerance,
            output_path             = NA,
            make_plots              = FALSE,
            verbose                 = FALSE,
            ...)
    ), silent = TRUE)
    if (inherits(res, "try-error")) {
      return(data.frame(n_plateaus = NA_integer_, ages = NA_character_,
                        mean_age = NA_real_))
    }
    ages <- res$summary[["Final age (Ma)"]]
    data.frame(
      n_plateaus = length(ages),
      ages       = if (length(ages)) paste(signif(ages, 6), collapse = "; ") else "",
      mean_age   = if (length(ages)) mean(ages) else NA_real_
    )
  }

  if (n_cores > 1L && requireNamespace("parallel", quietly = TRUE)) {
    cl <- parallel::makeCluster(min(n_cores, n))
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterExport(cl, c("adept", "ablate", "grid", "file_path",
                                  "calibration_uncertainty", "direction_method",
                                  "direction_tolerance"),
                            envir = environment())
    parallel::clusterEvalQ(cl, suppressMessages(library(ADEPT)))
    parts <- parallel::parLapply(cl, seq_len(n), ablate)
  } else {
    parts <- vector("list", n)
    for (k in seq_len(n)) {
      parts[[k]] <- ablate(k)
      if (is.function(progress)) {
        try(progress(k / n, sprintf("%d/%d", k, n)), silent = TRUE)
      }
      if (isTRUE(verbose)) {
        message(sprintf("[%d/%d] vt=%.4f, min_res=%s, %s -> %s plateau(s)",
                        k, n, grid$variance_threshold[k],
                        grid$min_plateau_resolution[k],
                        grid$filter_direction[k],
                        parts[[k]]$n_plateaus))
      }
    }
  }

  out <- cbind(grid, do.call(rbind, parts))
  row.names(out) <- NULL
  out
}

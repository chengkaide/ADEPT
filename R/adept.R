# ADEPT: Automated Depth Profiling Technique
# adept.R — Main exported function
#
# Wang et al. (2026) method for automated identification of age plateaus
# from LA-ICP-MS depth profiling data.
#
# This file holds the output schema and the entry point. The per-zircon work
# lives in zircon.R (grouping, ablation window, plateau fitting); the
# statistics live in segstats.R and processing.R; the filter cascade in
# filtering.R. `adept()` itself does three things only: read the input, walk
# over the zircons, and assemble the output.
#
# Dependencies: base R + (optional) writexl / mcp / shiny.
# Everything else — OOXML I/O, U-Pb age conversion, PELT segmentation,
# ARIMA order selection and plotting — is implemented inside the package.

# ---------------------------------------------------------------------------
# Output schema
# ---------------------------------------------------------------------------

ADEPT_BASE_COLS <- c(
  "Analysis", "Group", "Points",
  "Start", "End",
  "Integration time", "Max integration time", "Min integration time",
  "Standardized mean", "Segmentation mean",
  "Variance",
  "Calibration uncertainty", "Plateau uncertainty", "Total uncertainty",
  "Filter 1", "Filter 2", "Filter 3", "Filter 4",
  "Total integration time numbers", "Total segments", "Plateau serial numbers",
  "Final integration time numbers", "Final serial number",
  "Final age (Ma)", "Final total uncertainty (Ma)",
  "Final total uncertainty incl. decay (Ma)",
  "Random uncertainty (Ma)", "Systematic uncertainty (Ma)",
  "Decay constant uncertainty (Ma)", "Decay system",
  "Relative uncertainty (%)",
  "MSWD", "MSWD probability", "Points with uncertainty",
  "Confirmed plateaus",
  "Concordance (%)",
  "Pb206/U238 age mean (Ma)", "Pb206/U238 total uncertainty (Ma)",
  "Pb207/U235 age mean (Ma)", "Pb207/U235 total uncertainty (Ma)",
  "Pb207/Pb206 age mean (Ma)", "Pb207/Pb206 total uncertainty (Ma)"
)

ADEPT_MCMC_OUT_COLS <- c(
  "Analysis", "Group", "Points",
  "Start", "End",
  "Integration time", "Max integration time", "Min integration time",
  "Standardized mean", "Segmentation mean",
  "Standardized slope", "Loess slope",
  "Variance",
  "Calibration uncertainty", "Plateau uncertainty", "Total uncertainty",
  "Plateau serial numbers",
  "Filter 1", "Filter 2", "Filter 3", "Filter 4",
  "Final serial number", "Total integration time Numbers", "Total segments",
  "MCMC mean", "MCMC lower", "MCMC upper", "MCMC sigma",
  "Rhat", "MCMC n.eff", "MCMC integration time",
  "Filter MCMC", "Final integration time numbers",
  "Final age (Ma)", "Final total uncertainty (Ma)",
  "Final total uncertainty incl. decay (Ma)",
  "Random uncertainty (Ma)", "Systematic uncertainty (Ma)",
  "Decay constant uncertainty (Ma)", "Decay system",
  "Relative uncertainty (%)",
  "Confirmed plateaus",
  "Concordance (%)",
  "Pb206/U238 age mean (Ma)", "Pb206/U238 total uncertainty (Ma)",
  "Pb207/U235 age mean (Ma)", "Pb207/U235 total uncertainty (Ma)",
  "Pb207/Pb206 age mean (Ma)", "Pb207/Pb206 total uncertainty (Ma)"
)

ADEPT_SUMMARY_COLS <- c(
  "Analysis", "Group", "Points",
  "Final serial number",
  "Integration time",
  "Final age (Ma)", "Final total uncertainty (Ma)",
  "Final total uncertainty incl. decay (Ma)",
  "Random uncertainty (Ma)", "Systematic uncertainty (Ma)",
  "Decay constant uncertainty (Ma)", "Decay system",
  "Relative uncertainty (%)",
  "MSWD", "MSWD probability", "Points with uncertainty",
  "Confirmed plateaus",
  "Total segments",
  "Concordance (%)",
  "Pb206/U238 age mean (Ma)", "Pb206/U238 total uncertainty (Ma)",
  "Pb207/U235 age mean (Ma)", "Pb207/U235 total uncertainty (Ma)",
  "Pb207/Pb206 age mean (Ma)", "Pb207/Pb206 total uncertainty (Ma)"
)

#' Build one zircon's output block (vectorised over plateaus)
#'
#' @param seg Plateau table for one zircon.
#' @param analysis_name,group_counter,n_points Scalars replicated per row.
#' @param out_cols,target Column schema and the target column names.
#' @return data.frame with length(out_cols) columns and nrow(seg) rows.
#' @keywords internal
adept_block <- function(seg, analysis_name, group_counter, n_points,
                        out_cols, extra_mean_cols, mcmc) {
  m <- nrow(seg)
  col_of <- function(name) {
    if (is.null(name) || !(name %in% colnames(seg))) return(rep(NA_real_, m))
    v <- seg[[name]]
    if (is.null(v)) rep(NA_real_, m) else v
  }
  chr_of <- function(name) {
    if (is.null(name) || !(name %in% colnames(seg))) return(rep(NA_character_, m))
    as.character(seg[[name]])
  }

  vals <- vector("list", length(out_cols))
  names(vals) <- out_cols

  set <- function(col, v) vals[[col]] <<- v

  set("Analysis", rep(as.character(analysis_name), m))
  set("Group",    rep(as.numeric(group_counter), m))
  set("Points",   rep(as.numeric(n_points), m))
  set("Start",    col_of("Start"))
  set("End",      col_of("End"))
  set("Integration time",     col_of("Time_step"))
  set("Max integration time", col_of("Max_step"))
  set("Min integration time", col_of("Min_step"))
  set("Standardized mean",    col_of("standardized_Mean"))
  set("Segmentation mean",    col_of("Segment_Mean"))
  set("Variance",             col_of("Variance"))
  set("Calibration uncertainty", col_of("Calibration_uncertainty"))
  set("Plateau uncertainty",     col_of("Plateau_uncertainty"))
  set("Total uncertainty",       col_of("Total_uncertainty"))
  # v1.3.0: present only when the input carried per-point sigmas. col_of()
  # returns NA for an absent column, so the output shape does not change.
  set("MSWD",                    col_of("MSWD"))
  set("MSWD probability",        col_of("MSWD_prob"))
  set("Points with uncertainty", col_of("Uncertainty_n"))
  set("Filter 1", col_of("Filter_1"))
  set("Filter 2", col_of("Filter_2"))
  set("Filter 3", col_of("Filter_3"))
  set("Filter 4", col_of("Filter_4"))

  if (isTRUE(mcmc)) {
    set("Standardized slope", col_of("standardized_slope"))
    set("Loess slope",        col_of("loess_slope"))
    set("Plateau serial numbers", col_of("Number"))
    set("Final serial number",    col_of("Final_Serial_Number"))
    set("Total integration time Numbers", col_of("Plateau_Numbers"))
    set("MCMC mean",  col_of("mean"))
    set("MCMC lower", col_of("lower"))
    set("MCMC upper", col_of("upper"))
    set("MCMC sigma", col_of("sigma"))
    set("Rhat",       col_of("Rhat"))
    set("MCMC n.eff", col_of("n.eff"))
    set("MCMC integration time", col_of("mcp_step"))
    set("Filter MCMC", col_of("Filter_mcp"))
    set("Final integration time numbers", col_of("Final_step_Numbers"))
  } else {
    set("Total integration time numbers", col_of("Plateau_Numbers"))
    set("Plateau serial numbers",         col_of("Number"))
    set("Final integration time numbers", col_of("Final_step_Numbers"))
    set("Final serial number",            col_of("Final_Serial_Number"))
  }

  set("Final age (Ma)",               col_of("Final_Age"))
  set("Final total uncertainty (Ma)", col_of("Final_total_uncertainty"))
  set("Final total uncertainty incl. decay (Ma)",
      col_of("Final_total_uncertainty_full"))
  set("Random uncertainty (Ma)",          col_of("Random_uncertainty"))
  set("Systematic uncertainty (Ma)",      col_of("Systematic_uncertainty"))
  set("Decay constant uncertainty (Ma)",  col_of("Decay_constant_uncertainty"))
  set("Decay system",                     chr_of("Decay_system"))
  set("Relative uncertainty (%)",         col_of("Relative_uncertainty_pct"))
  set("Confirmed plateaus",               col_of("Confirmed_Plateaus"))
  set("Total segments",                   col_of("Total_Segments"))
  set("Concordance (%)",              col_of("Concordance"))
  set("Pb206/U238 age mean (Ma)",           col_of("Age68_Mean"))
  set("Pb206/U238 total uncertainty (Ma)",  col_of("Age68_Total_uncertainty"))
  set("Pb207/U235 age mean (Ma)",           col_of("Age75_Mean"))
  set("Pb207/U235 total uncertainty (Ma)",  col_of("Age75_Total_uncertainty"))
  set("Pb207/Pb206 age mean (Ma)",          col_of("Age76_Mean"))
  set("Pb207/Pb206 total uncertainty (Ma)", col_of("Age76_Total_uncertainty"))

  for (nm in extra_mean_cols) {
    vals[[nm]] <- col_of(nm)
  }

  as.data.frame(vals, stringsAsFactors = FALSE, optional = TRUE)
}

#' Empty output block (one row of NA) for a zircon that cannot be processed
#' @keywords internal
adept_blank_block <- function(analysis_name, group_counter, n_points, out_cols) {
  vals <- as.list(rep(NA_real_, length(out_cols)))
  names(vals) <- out_cols
  vals[["Analysis"]] <- as.character(analysis_name)
  vals[["Group"]]    <- as.numeric(group_counter)
  vals[["Points"]]   <- as.numeric(n_points)
  as.data.frame(vals, stringsAsFactors = FALSE, optional = TRUE)
}

#' Automated Depth Profiling Technique
#'
#' Quantitatively identifies and extracts age plateaus from LA-ICP-MS
#' U-Pb depth profiling data. Supports four input formats: direct ages
#' (Age68/Age75/Age76), raw isotope counts (Pb206/Pb207/U238), isotopic ratios
#' (Pb206_U238/Pb207_U235/Pb207_Pb206), or already-reduced ages with per-point
#' uncertainties (Age68 plus Age68_1s).
#'
#' Each zircon is processed independently: the rows sharing an Analysis value
#' form one unit, the effective-ablation-time window is applied, the profile is
#' screened and smoothed, PELT finds the segment boundaries, and a four-step
#' cascade keeps the plateaus that pass. When the input carries per-point
#' 1-sigmas the plateau age becomes an inverse-variance weighted mean and MSWD
#' is reported alongside it.
#'
#' @param file_path Path to the input Excel (.xlsx) or CSV file.
#' @param chunk_size Number of rows per processing chunk. Only used for input
#'   with no Analysis column, where the row count is the only grouping signal;
#'   otherwise consecutive rows sharing an Analysis value form one zircon.
#' @param lower_ablation_time Minimum effective ablation time in seconds (default 29)
#' @param upper_ablation_time Maximum effective ablation time in seconds (default 58)
#' @param max_age_limit Maximum valid age in Ma (default 4540)
#' @param min_age_limit Minimum valid age in Ma (default 0)
#' @param min_plateau_resolution Minimum plateau duration in seconds.
#'   \code{NULL} or values < 5 default to 5 seconds.
#' @param variance_threshold Maximum allowed intra-plateau variance (default 0.1192)
#' @param filter_direction \code{"Forward"} keeps ascending age sequences;
#'   \code{"Reverse"} keeps descending age sequences.
#' @param direction_method,direction_tolerance How the directional step selects
#'   plateaus; see \code{\link{filter_direction}}.
#' @param outlier_method Outlier detection before smoothing:
#'   \code{"arima"} (default, matches the published method), \code{"mad"}
#'   (fast robust z-score) or \code{"none"}.
#' @param outlier_sd Threshold in residual standard deviations (default 2).
#' @param preprocess \code{"arima_loess"} (default, the published pipeline) or
#'   \code{"robust_loess"}, which skips the ARIMA screen and fits a robust
#'   M-estimator LOESS instead.
#' @param smooth \code{"loess"} (default) or \code{"none"}. Use \code{"none"}
#'   for input that already carries a down-hole fractionation correction: such
#'   a profile is flat inside a domain, and LOESS would round off the real
#'   domain boundaries.
#' @param calibration_uncertainty Relative 1-sigma reproducibility of the
#'   primary reference material (default 0.03, i.e. 3 \%). Set it to \code{0}
#'   when the input's own per-point sigmas already include the external
#'   reproducibility, or the term is counted twice.
#' @param u238u235 238U/235U ratio used for count-based input (default 137.818).
#' @param validate_input Logical. Report structural input problems up front?
#' @param mcmc Logical. Run Bayesian MCMC posterior analysis? (default \code{FALSE}).
#'   Requires the suggested package \pkg{mcp} (and JAGS); if unavailable the
#'   MCMC columns are \code{NA} rather than an error.
#' @param make_plots Logical. Build depth profile plots? (default \code{TRUE})
#' @param save_plots_to_disk Logical. Write plot PDFs to \code{plot_dir}?
#' @param output_path Output Excel path. \code{NULL} auto-generates from the
#'   input name; \code{NA} skips writing the workbook.
#' @param plot_dir Directory for plot PDFs. \code{NULL} uses the input directory.
#' @param keep_profiles Logical. Also return per-zircon series and plateau
#'   tables under \code{$profiles}? (default \code{FALSE})
#' @param progress Optional callback \code{function(fraction, detail)}.
#' @param verbose Logical. Emit per-zircon messages? (default \code{TRUE})
#' @param plot Deprecated alias for \code{make_plots}.
#'
#' @return Invisibly returns a list with \code{summary}, \code{full},
#'   \code{plots} and (optionally) \code{profiles}.
#'
#' @seealso \code{\link{adept_sensitivity}} to quantify how much the result
#'   depends on \code{variance_threshold} and \code{min_plateau_resolution}.
#'
#' @export
#' @importFrom stats loess predict var sd residuals lm coef complete.cases median mad setNames
#' @importFrom utils head
#'
#' @examples
#' \dontrun{
#' result <- adept("Input.xlsx")
#' result <- adept("Input.xlsx", filter_direction = "Reverse")
#' result <- adept("Input.xlsx", mcmc = TRUE)
#' }
adept <- function(
    file_path,
    chunk_size                = 411,
    lower_ablation_time       = 29,
    upper_ablation_time       = 58,
    max_age_limit             = 4540,
    min_age_limit             = 0,
    min_plateau_resolution    = NULL,
    variance_threshold        = 0.1192,
    filter_direction          = c("Forward", "Reverse"),
    direction_method          = c("monotonic", "strict"),
    direction_tolerance       = 0.02,
    outlier_method            = c("arima", "mad", "none"),
    outlier_sd                = 2,
    preprocess                = c("arima_loess", "robust_loess"),
    smooth                    = c("loess", "none"),
    calibration_uncertainty   = 0.03,
    u238u235                  = ADEPT_U238U235,
    validate_input            = TRUE,
    mcmc                      = FALSE,
    make_plots                = TRUE,
    save_plots_to_disk        = TRUE,
    output_path               = NULL,
    plot_dir                  = NULL,
    keep_profiles             = FALSE,
    progress                  = NULL,
    verbose                   = TRUE,
    plot                      = TRUE
) {
  filter_direction  <- match.arg(filter_direction)
  direction_method  <- match.arg(direction_method)
  outlier_method    <- match.arg(outlier_method)
  preprocess        <- match.arg(preprocess)
  smooth            <- match.arg(smooth)

  if (!is.numeric(calibration_uncertainty) ||
      length(calibration_uncertainty) != 1L ||
      !is.finite(calibration_uncertainty) || calibration_uncertainty < 0) {
    stop("`calibration_uncertainty` must be a single non-negative number ",
         "(it is a relative 1-sigma, e.g. 0.03 for 3 %).", call. = FALSE)
  }

  if (!missing(plot) && missing(make_plots)) make_plots <- isTRUE(plot)
  if (!is.null(progress) && !is.function(progress)) {
    stop("`progress` must be a function or NULL.", call. = FALSE)
  }
  report <- function(fraction, detail) {
    if (is.function(progress)) try(progress(fraction, detail), silent = TRUE)
    # v1.3.0: also surface it under verbose. Without this a zircon that drops
    # out early ("too few points", "empty window", ...) disappears silently,
    # and the only symptom is a short summary table.
    if (isTRUE(verbose)) message("  ", detail)
    invisible(NULL)
  }

  if (!file.exists(file_path)) {
    stop("Input file not found: ", file_path, call. = FALSE)
  }
  if (!is.numeric(u238u235) || length(u238u235) != 1L ||
      !is.finite(u238u235) || u238u235 <= 0) {
    stop("`u238u235` must be a single positive number.", call. = FALSE)
  }

  data_list   <- read_input(file_path)
  sheet_names <- names(data_list)

  blocks    <- list()
  plots     <- list()
  profiles  <- list()
  out_cols  <- if (isTRUE(mcmc)) ADEPT_MCMC_OUT_COLS else ADEPT_BASE_COLS

  # ---- Pre-scan: master list of extra numeric columns ----------------------
  all_extra_names <- character(0)
  for (sheet_name in sheet_names) {
    raw_sheet <- data_list[[sheet_name]]
    all_extra_names <- union(
      all_extra_names,
      detect_extra_names(raw_sheet, 1, min(20L, nrow(raw_sheet)))
    )
  }
  extra_mean_cols <- if (length(all_extra_names) > 0) {
    paste0(all_extra_names, "_Mean")
  } else character(0)
  out_cols <- c(out_cols, extra_mean_cols)

  # ---- Output path ---------------------------------------------------------
  base <- sub("\\.[^.]*$", "", basename(file_path))
  if (is.null(output_path)) {
    output_path <- file.path(dirname(file_path), paste0(base, "_Output.xlsx"))
  }
  if (is.null(plot_dir)) plot_dir <- dirname(file_path)

  # ---- Grouping and progress bookkeeping -----------------------------------
  # One unit of work is one zircon. The row ranges are computed once and reused
  # by the main loop, so the progress denominator and the loop cannot disagree
  # about how many zircons there are.
  groups <- setNames(lapply(sheet_names, function(s) {
    adept_zircon_groups(data_list[[s]], chunk_size, s)
  }), sheet_names)
  total_units <- sum(vapply(groups, function(g) as.numeric(g$n), numeric(1)))
  if (total_units < 1) total_units <- 1
  done_units <- 0

  # Everything the per-zircon steps need, in one object. Each step then takes
  # (data, cfg) instead of a dozen positional arguments, and adding a parameter
  # to adept() means adding one line here rather than threading it through
  # four call signatures.
  cfg <- list(
    u238u235                = u238u235,
    validate_input          = validate_input,
    lower_ablation_time     = lower_ablation_time,
    upper_ablation_time     = upper_ablation_time,
    preprocess              = preprocess,
    smooth                  = smooth,
    outlier_method          = outlier_method,
    outlier_sd              = outlier_sd,
    calibration_uncertainty = calibration_uncertainty,
    min_age_limit           = min_age_limit,
    max_age_limit           = max_age_limit,
    variance_threshold      = variance_threshold,
    min_plateau_resolution  = min_plateau_resolution,
    filter_direction        = filter_direction,
    direction_method        = direction_method,
    direction_tolerance     = direction_tolerance,
    all_extra_names         = all_extra_names
  )

  # ---- Main loop (one iteration per zircon) --------------------------------
  for (sheet_name in sheet_names) {
    raw <- data_list[[sheet_name]]
    grp <- groups[[sheet_name]]

    if (isTRUE(verbose)) {
      message(sprintf("Processing: '%s' (%d zircon(s))", sheet_name, grp$n))
    }

    group_counter <- 1
    for (i in seq_len(grp$n)) {
      start_row <- grp$starts[i]
      end_row   <- grp$ends[i]

      r <- adept_one_zircon(raw, start_row, end_row, sheet_name, cfg)

      if (!isTRUE(r$ok)) {
        # A zircon that cannot be processed still occupies a row in the output,
        # so that the summary accounts for every unit the progress bar counted.
        blocks[[length(blocks) + 1L]] <-
          adept_blank_block(grp$labels[i], group_counter, r$points, out_cols)
        done_units <- done_units + 1
        report(done_units / total_units,
               sprintf("%s #%d (%s)", sheet_name, i, r$label))
        group_counter <- group_counter + 1
        next
      }

      n_confirmed <- r$n_confirmed
      if (isTRUE(verbose)) {
        message(sprintf("[%d/%d] %s - %d plateau(s) confirmed (%s)",
                        i, grp$n, r$analysis, n_confirmed, filter_direction))
      }

      segments <- r$segments
      if (isTRUE(mcmc)) segments <- run_mcmc(r$data, segments)

      if (isTRUE(make_plots)) {
        plots[[length(plots) + 1L]] <- plot_depth_profile(
          r$data, segments,
          title = paste("Analysis:", r$analysis),
          label = paste0(group_counter, "_", r$analysis)
        )
      }

      if (isTRUE(keep_profiles)) {
        profiles[[length(profiles) + 1L]] <- list(
          sheet    = sheet_name,
          group    = group_counter,
          analysis = r$analysis,
          data     = r$data,
          segments = segments
        )
      }

      blocks[[length(blocks) + 1L]] <- adept_block(
        segments,
        analysis_name = r$analysis,
        group_counter = group_counter,
        n_points      = nrow(r$data),
        out_cols      = out_cols,
        extra_mean_cols = extra_mean_cols,
        mcmc          = mcmc
      )

      done_units <- done_units + 1
      report(done_units / total_units,
             sprintf("%s #%d - %d plateau(s)", sheet_name, i, n_confirmed))
      group_counter <- group_counter + 1
    }
  }

  report(1, "Finalising output")

  # ---- Assemble -------------------------------------------------------------
  if (length(blocks) > 0) {
    output <- do.call(rbind, blocks)
    row.names(output) <- NULL
    colnames(output)  <- out_cols
  } else {
    warning("No plateaus were produced. The input may be empty or every ",
            "zircon failed the quality gates.", call. = FALSE)
    output <- as.data.frame(
      setNames(lapply(out_cols, function(cn)
        if (cn == "Analysis") character(0) else numeric(0)), out_cols),
      stringsAsFactors = FALSE, optional = TRUE)
  }

  simplified_cols <- intersect(c(ADEPT_SUMMARY_COLS, extra_mean_cols),
                               colnames(output))
  output_simple <- output[, simplified_cols, drop = FALSE]
  if (nrow(output_simple) > 0) {
    output_simple <- output_simple[!is.na(output_simple[["Final age (Ma)"]]), ,
                                   drop = FALSE]
    row.names(output_simple) <- NULL
  }

  if (!is.null(output_path) && !is.na(output_path)) {
    write_xlsx_base(list(Summary = output_simple, Full_Results = output),
                    output_path)
    if (isTRUE(verbose)) message("Results written to: ", output_path)
  }

  if (isTRUE(make_plots) && isTRUE(save_plots_to_disk) && length(plots) > 0) {
    save_plots(plots, plot_dir)
  }

  res <- list(summary = output_simple, full = output, plots = plots)
  if (isTRUE(keep_profiles)) res$profiles <- profiles
  invisible(res)
}

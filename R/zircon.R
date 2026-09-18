# ADEPT: Automated Depth Profiling Technique
# zircon.R - Per-zircon pipeline steps
#
# `adept()` is the only entry point, but it is not where the work happens.
# Everything below is one named step of the per-zircon pipeline, so that a
# change to (say) the ablation window cannot silently also change the
# filtering. Each step is a pure function: it takes data plus a configuration
# list, and returns data plus a status. That is what makes the steps testable
# on their own, and it is what stops the same check from being written twice
# in two different places with two different behaviours.
#
# Deliberately NOT here: bookkeeping. Which blocks, plots and profiles have
# been produced so far, and where the progress bar is, belongs to `adept()`
# alone. Each step below reports failure as `ok = FALSE` plus a short label,
# and lets the caller decide what to record.
#
# One zircon goes through four stages, in this order:
#
#   adept_zircon_groups()   split a sheet into per-zircon row ranges
#   adept_one_zircon()      parse -> prepare -> fit, one zircon
#     adept_prepare_window()  pick the ablation window, build Raw_Age (+ sigma)
#     adept_fit_plateaus()    screen -> smooth -> PELT -> statistics -> filter

#' Split a sheet into per-zircon row ranges
#'
#' The unit of processing is one zircon, and the unit of output is one row per
#' plateau of one zircon. The row ranges therefore have to follow the Analysis
#' column rather than the row count: a sheet that stacks several zircons must
#' not be cut mid-zircon, and two zircons must not be merged into one block.
#' Consecutive rows sharing an Analysis value form one unit.
#'
#' `chunk_size` is only the fallback for input with no Analysis column at all,
#' where the row count is the only grouping signal available.
#'
#' @param d One sheet, as read from the input file.
#' @param chunk_size Rows per group when there is no Analysis column.
#' @param sheet_name Sheet name, used to build labels in the fallback case.
#' @return List with `starts`, `ends`, `labels`, `n` (number of zircons) and
#'   `analysis_col` (`NA` when the sheet has no Analysis column).
#' @keywords internal
adept_zircon_groups <- function(d, chunk_size, sheet_name) {
  an_col <- intersect(c("Analysis", "Analysis_"), colnames(d))[1]

  if (!is.na(an_col)) {
    runs   <- rle(as.character(d[[an_col]]))
    ends   <- cumsum(runs$lengths)
    # An empty sheet has no runs, which would leave `ends` empty while `starts`
    # still held one element. One group covering no rows is the honest reading.
    if (length(ends) == 0L) ends <- 0L
    starts <- c(1L, utils::head(ends, -1L) + 1L)
    labels <- as.character(d[[an_col]][starts])
  } else {
    k      <- max(1L, ceiling(nrow(d) / chunk_size))
    starts <- (seq_len(k) - 1L) * chunk_size + 1L
    ends   <- pmin(seq_len(k) * chunk_size, nrow(d))
    labels <- paste0(sheet_name, "_", starts)
  }

  list(starts = starts, ends = ends, labels = labels,
       n = length(starts), analysis_col = an_col)
}

#' Does this table still contain a usable age?
#'
#' The 207Pb/235U and 207Pb/206Pb ages are optional (a Format 4 input from an
#' external reduction carries neither), so only the columns that are present
#' are tested. Testing an absent column would hand `all(is.na(NULL))` — which
#' is TRUE — straight into the guard, and every zircon would be rejected as
#' having no ages. That mistake is what made the v1.2.0 code unable to read a
#' two-column input at all.
#'
#' @param d data.frame with an Age68 column.
#' @param has_75,has_76 Whether the optional age columns are present.
#' @return TRUE when at least one usable age column remains.
#' @keywords internal
adept_has_ages <- function(d, has_75, has_76) {
  !(all(is.na(d$Age68)) ||
    (has_75 && all(is.na(d$Age75))) ||
    (has_76 && all(is.na(d$Age76))))
}

#' Select the ablation window and build the age series used downstream
#'
#' Applies the effective-ablation-time window, then builds `Raw_Age` — the
#' series the plateau statistics are computed on — plus `Raw_Age_sigma`
#' alongside it when the input carried per-point 1-sigmas.
#'
#' The window is applied to every input, including a pre-reduced one. An
#' external reduction that has already cut its own window should pass
#' `lower_ablation_time = 0` and an upper bound past its last window, rather
#' than have ADEPT guess: a per-point sigma column says nothing about where the
#' window ends, and guessing wrong silently changes how many points enter the
#' segmentation.
#'
#' @param segment.data Parsed table for one zircon, with `.ROWID.` attached.
#' @param cfg Pipeline configuration (see `adept()`).
#' @param sheet_name,start_row,end_row Used only to label warnings.
#' @return List: on success `ok = TRUE` with `data`, `has_75`, `has_76`; on
#'   failure `ok = FALSE` with `points` and a short `label` naming the reason.
#' @keywords internal
adept_prepare_window <- function(segment.data, cfg, sheet_name, start_row,
                                end_row) {
  fail <- function(points, label) list(ok = FALSE, points = points, label = label)

  has_75 <- "Age75" %in% names(segment.data)
  has_76 <- "Age76" %in% names(segment.data)

  if (!adept_has_ages(segment.data, has_75, has_76)) {
    return(fail(0, "no ages"))
  }

  num_cols <- setdiff(names(segment.data), c("Analysis", ".ROWID."))
  for (cn in num_cols) {
    segment.data[[cn]] <- suppressWarnings(as.numeric(segment.data[[cn]]))
  }
  need <- c("Analysis", "Time", "Age68",
            if (has_75) "Age75", if (has_76) "Age76")
  keep <- complete.cases(segment.data[, need, drop = FALSE])
  segment.data <- segment.data[keep, , drop = FALSE]

  in_win      <- which(segment.data$Time >= cfg$lower_ablation_time &
                       segment.data$Time <= cfg$upper_ablation_time)
  subset_data <- segment.data[in_win, , drop = FALSE]

  if (nrow(subset_data) == 0) {
    return(fail(0, "empty window"))
  }

  if (isTRUE(cfg$validate_input)) {
    problems <- validate_segment_data(
      subset_data,
      label = sprintf("'%s' rows %d-%d (ablation window)", sheet_name,
                      start_row, end_row),
      check_values = TRUE)
    # Non-fatal: deep-profile files routinely contain nonsense ages in the
    # gas-blank and wash-out portions, and those rows are already excluded.
    if (length(problems) > 0L) {
      warning(sprintf("Sheet '%s' rows %d-%d (ablation window): %s",
                      sheet_name, start_row, end_row,
                      paste(problems, collapse = "; ")), call. = FALSE)
    }
  }

  # 207Pb/206Pb is optional, so use a safe stand-in for the branch test
  # instead of letting ifelse() receive a zero-length vector.
  a76 <- if (is.null(subset_data$Age76)) {
    rep(NA_real_, nrow(subset_data))
  } else {
    subset_data$Age76
  }

  # Invariant from here on: complete.cases() below drops every row with an NA
  # in a required age column, so Age68 cannot be all-NA inside the window, and
  # a sheet whose Age68 is entirely blank never reaches this function at all
  # (format detection refuses it). A second "no ages" guard after the window
  # was therefore dead code, and v1.3.0 briefly carried two copies of it that
  # could drift apart. There is one guard now, and it is the one above.
  subset_data$Raw_Age <- ifelse(subset_data$Age68 < 1000,
                                subset_data$Age68,
                                ifelse(a76 > 1000, a76, NA))

  # The sigma branches mirror the Raw_Age branches exactly, so the sigma always
  # belongs to whichever age Raw_Age was taken from.
  if (!is.null(subset_data$Age68_1s)) {
    s68 <- subset_data$Age68_1s
    s76 <- if (!is.null(subset_data$Age76_1s)) {
      subset_data$Age76_1s
    } else {
      rep(NA_real_, nrow(subset_data))
    }
    subset_data$Raw_Age_sigma <-
      ifelse(subset_data$Age68 < 1000, s68,
             ifelse(a76 > 1000, s76, NA))
  }

  list(ok = TRUE, data = subset_data, has_75 = has_75, has_76 = has_76)
}

#' Screen, smooth, segment and filter one zircon
#'
#' The whole statistical chain, from the ablation window to the four-step
#' filter cascade, in the order the method defines it. Everything here is
#' vectorised over segments; there is no per-point R loop.
#'
#' `preprocess = "arima_loess"` is the published pipeline: an ARIMA residual
#' screen, then an ordinary LOESS fit. `"robust_loess"` skips the screen and
#' lets a robust M-estimator LOESS downweight outliers instead — fewer steps,
#' no model-order selection, no hard deletion of points.
#'
#' `smooth = "none"` is for input that already carries a down-hole
#' fractionation correction: the profile is then flat inside a domain, and
#' LOESS would round off the real domain boundaries.
#'
#' @param subset_data Output of `adept_prepare_window()`.
#' @param extra_raw,extra_names Extra (non-age) columns carried alongside.
#' @param cfg Pipeline configuration (see `adept()`).
#' @param has_75,has_76 Whether the optional age columns are present.
#' @return List: on success `ok = TRUE` with `data`, `segments`,
#'   `extra_names` and `n_confirmed`; on failure `ok = FALSE` with `points`
#'   and a short `label` naming the reason.
#' @keywords internal
adept_fit_plateaus <- function(subset_data, extra_raw, extra_names, cfg,
                               has_75, has_76) {
  fail <- function(points, label) list(ok = FALSE, points = points, label = label)

  if (identical(cfg$preprocess, "arima_loess")) {
    subset_data$Age68 <- arima_outlier(subset_data$Age68, cfg$outlier_method,
                                       cfg$outlier_sd)
    # Only screen the optional columns that are actually present.
    if (has_75) {
      subset_data$Age75 <- arima_outlier(subset_data$Age75, cfg$outlier_method,
                                         cfg$outlier_sd)
    }
    if (has_76) {
      subset_data$Age76 <- arima_outlier(subset_data$Age76, cfg$outlier_method,
                                         cfg$outlier_sd)
    }
  }

  subset_data <- discordance_filter(subset_data)
  subset_data$subset_Age68 <- mean_fill(subset_data$subset_Age68,
                                        subset_data$Age68)
  if (has_76) {
    subset_data$subset_Age76 <- mean_fill(subset_data$subset_Age76,
                                          subset_data$Age76)
  } else {
    subset_data$subset_Age76 <- rep(NA_real_, nrow(subset_data))
  }
  subset_data$subset_Age <- ifelse(subset_data$subset_Age68 < 1000,
                                   subset_data$subset_Age68,
                                   ifelse(subset_data$subset_Age76 > 1000,
                                          subset_data$subset_Age76, NA))

  ok_age <- !is.na(subset_data$subset_Age)
  if (sum(ok_age) < 10) {
    return(fail(nrow(subset_data), "too few points"))
  }
  subset_data <- subset_data[ok_age, , drop = FALSE]
  subset_data$Row_Number <- seq_len(nrow(subset_data))

  # ---- Merge extra columns -------------------------------------------------
  if (length(extra_names) > 0 && !is.null(extra_raw)) {
    matched_rows <- match(subset_data$.ROWID., extra_raw$.ROWID.)
    extra_names  <- intersect(extra_names, cfg$all_extra_names)
    for (col in extra_names) {
      subset_data[[col]] <- suppressWarnings(
        as.numeric(as.character(extra_raw[[col]][matched_rows])))
    }
  } else {
    extra_names <- character(0)
  }

  # ---- Smoothing -----------------------------------------------------------
  subset_data <- loess_segment(
    subset_data, span = 0.15,
    family = if (identical(cfg$preprocess, "robust_loess")) "symmetric"
             else "gaussian",
    smooth = cfg$smooth)
  if (nrow(subset_data) == 0 ||
      sum(!is.na(subset_data$standardized_loess)) < 10) {
    return(fail(nrow(subset_data), "smoothing failed"))
  }

  # ---- PELT ----------------------------------------------------------------
  pelt         <- pelt_segmentation(subset_data$standardized_loess,
                                    nrow(subset_data))
  changepoints <- pelt$changepoints
  subset_data$Change_loess <- ifelse(subset_data$Row_Number %in% changepoints,
                                     "YES", "NO")

  seg_result  <- build_segments(subset_data, changepoints)
  segments    <- seg_result$segments
  subset_data <- seg_result$df
  seg_starts  <- seg_result$segment_starts
  seg_ends    <- seg_result$segment_ends

  if (nrow(segments) == 0) {
    return(fail(nrow(subset_data), "no segments"))
  }

  # ---- Plateau statistics (vectorised) -------------------------------------
  segments <- calc_slopes(subset_data, segments, seg_starts, seg_ends)
  segments$Variance <- calc_variance(subset_data, seg_starts, seg_ends)
  segments <- calc_uncertainty(segments, subset_data, seg_starts, seg_ends,
                               calibration_uncertainty =
                                 cfg$calibration_uncertainty)
  segments <- calc_extra_means(segments, subset_data, extra_names,
                               seg_starts, seg_ends)
  segments <- calc_age_means(segments, subset_data, seg_starts, seg_ends,
                             calibration_uncertainty =
                               cfg$calibration_uncertainty)
  segments <- calc_concordance(segments)

  segments <- apply_filters(segments,
                            min_age       = cfg$min_age_limit,
                            max_age       = cfg$max_age_limit,
                            var_threshold = cfg$variance_threshold,
                            min_res       = cfg$min_plateau_resolution,
                            direction     = cfg$filter_direction,
                            direction_method    = cfg$direction_method,
                            direction_tolerance = cfg$direction_tolerance)

  list(ok          = TRUE,
       data        = subset_data,
       segments    = segments,
       extra_names = extra_names,
       n_confirmed = sum(!is.na(segments$Filter_4)))
}

#' Run the full pipeline for one zircon
#'
#' Parse, prepare, fit. Returns a single uniform result so that `adept()` has
#' one branch to handle instead of eight: every way a zircon can drop out
#' (unparseable rows, failed validation, too few points, ...) comes back as
#' `ok = FALSE` with a short label, which becomes the progress message.
#'
#' @param raw The whole sheet, unparsed.
#' @param start_row,end_row Row range of this zircon.
#' @param sheet_name Sheet the rows came from, used to label warnings.
#' @param cfg Pipeline configuration (see `adept()`).
#' @return List with `ok`; when `ok` is FALSE also `points` and `label`; when
#'   TRUE also `data`, `segments`, `n_confirmed` and `analysis`.
#' @keywords internal
adept_one_zircon <- function(raw, start_row, end_row, sheet_name, cfg) {
  fail <- function(points, label) list(ok = FALSE, points = points, label = label)

  parsed <- try(parse_segment_data(raw, start_row, end_row, cfg$u238u235),
                silent = TRUE)
  if (inherits(parsed, "try-error")) {
    warning(sprintf("Sheet '%s' rows %d-%d skipped: %s", sheet_name,
                    start_row, end_row,
                    conditionMessage(attr(parsed, "condition"))),
            call. = FALSE)
    return(fail(0, "skipped"))
  }

  segment.data <- parsed$data
  extra_names  <- parsed$extra_names
  extra_raw    <- parsed$extra_data

  if (isTRUE(cfg$validate_input)) {
    problems <- validate_segment_data(
      segment.data,
      label = sprintf("'%s' rows %d-%d", sheet_name, start_row, end_row))
    # Only structural problems are fatal: a missing age column or a Time axis
    # that looks like milliseconds means the pipeline cannot run.
    fatal <- grepl("no recognisable|not monotonically|not numeric|is not numeric",
                   problems)
    if (any(fatal)) {
      warning(sprintf("Sheet '%s' rows %d-%d: %s", sheet_name, start_row,
                      end_row, paste(problems[fatal], collapse = "; ")),
              call. = FALSE)
      return(fail(0, "failed validation"))
    }
  }

  segment.data$.ROWID. <- seq_len(nrow(segment.data))
  if (length(extra_names) > 0 && !is.null(extra_raw)) {
    extra_raw$.ROWID. <- seq_len(nrow(extra_raw))
  }

  prep <- adept_prepare_window(segment.data, cfg, sheet_name, start_row,
                               end_row)
  if (!isTRUE(prep$ok)) return(fail(prep$points, prep$label))

  fit <- adept_fit_plateaus(prep$data, extra_raw, extra_names, cfg,
                            has_75 = prep$has_75, has_76 = prep$has_76)
  if (!isTRUE(fit$ok)) return(fail(fit$points, fit$label))

  list(ok          = TRUE,
       data        = fit$data,
       segments    = fit$segments,
       n_confirmed = fit$n_confirmed,
       analysis    = as.character(fit$data[1, "Analysis"]))
}

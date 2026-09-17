# ADEPT: Automated Depth Profiling Technique
# processing.R — Core data processing pipeline
#
# All statistics here are vectorised over segments (see segstats.R); the only
# remaining R-level loop in the package is the one that walks over zircons,
# which is inherent to the method (each zircon is fitted independently).
#
# External dependencies: none (base R only).

# ---------------------------------------------------------------------------
#  Input
# ---------------------------------------------------------------------------

#' Read all sheets from an Excel input file
#'
#' Uses the built-in OOXML reader (xlsx.R) so no Excel package is required.
#'
#' @param file_path Path to the input Excel file.
#' @return A named list of data.frames, one per sheet.
#' @keywords internal
read_input <- function(file_path) {
  if (!file.exists(file_path)) {
    stop("Input file not found: ", file_path, call. = FALSE)
  }
  if (grepl("\\.csv$", file_path, ignore.case = TRUE)) {
    return(list(Sheet1 = utils::read.csv(file_path, check.names = FALSE,
                                         stringsAsFactors = FALSE)))
  }
  if (grepl("\\.xls$", file_path, ignore.case = TRUE)) {
    stop("Legacy binary .xls files are not supported. ",
         "Please save the workbook as .xlsx or .csv.", call. = FALSE)
  }
  read_xlsx_base(file_path)
}

#' Detect "extra" numeric columns without converting ages
#'
#' Used by `adept()` for the pre-scan that builds the master list of extra
#' columns. It never computes ages, so it cannot fail on a sheet whose first
#' rows are blank.
#'
#' @param raw Raw data.frame from one sheet
#' @param start_row,end_row Row window to inspect
#' @return Character vector of extra numeric column names
#' @keywords internal
detect_extra_names <- function(raw, start_row, end_row) {
  colnames(raw) <- gsub(" |/", "_", colnames(raw))

  core_cols <- c("Analysis", "Time",
                 "Age68", "Age75", "Age76",
                 "Pb206", "Pb207", "U235", "U238",
                 "Pb206U238", "Pb207U235", "Pb207Pb206",
                 "Pb206_U238", "Pb207_U235", "Pb207_Pb206")

  cand <- setdiff(colnames(raw), core_cols)
  if (length(cand) == 0) return(character(0))

  end_row <- min(end_row, nrow(raw))
  if (end_row < start_row) return(character(0))
  sub <- raw[start_row:end_row, cand, drop = FALSE]

  ok <- vapply(cand, function(cn) {
    v <- sub[[cn]]
    if (is.numeric(v)) return(TRUE)
    cv <- suppressWarnings(as.numeric(as.character(v)))
    if (length(cv) == 0) return(FALSE)
    sum(!is.na(cv)) / length(cv) > 0.5
  }, logical(1), USE.NAMES = FALSE)

  cand[ok]
}

#' Parse raw segment data and extract age columns
#'
#' Detects the input format (Age68, Pb206, or Pb206_U238) and creates a
#' standardised data.frame with Analysis, Time, Age68, Age75, Age76.
#' Ages are derived with the built-in decay equations (isotopes.R).
#'
#' @param raw Raw data.frame from one sheet
#' @param start_row,end_row Row window for this chunk
#' @param u238u235 238U/235U ratio for count-based input
#' @return A list with `data`, `extra_names` and `extra_data`
#' @keywords internal
parse_segment_data <- function(raw, start_row, end_row, u238u235 = ADEPT_U238U235) {
  colnames(raw) <- gsub(" |/", "_", colnames(raw))

  if ("Age68" %in% colnames(raw) &&
      !all(is.na(raw$Age68[start_row:end_row]))) {
    segment.data <- data.frame(
      Analysis = raw$Analysis[start_row:end_row],
      Time     = raw$Time[start_row:end_row],
      Age68    = raw$Age68[start_row:end_row],
      Age75    = raw$Age75[start_row:end_row],
      Age76    = raw$Age76[start_row:end_row],
      stringsAsFactors = FALSE
    )
    segment.data[, 2:5] <- lapply(segment.data[, 2:5],
                                  function(v) suppressWarnings(as.numeric(v)))

  } else if ("Pb206" %in% colnames(raw) &&
             !all(is.na(raw$Pb206[start_row:end_row]))) {
    segment.data <- data.frame(
      Analysis = raw$Analysis[start_row:end_row],
      Time     = raw$Time[start_row:end_row],
      Pb206    = raw$Pb206[start_row:end_row],
      Pb207    = raw$Pb207[start_row:end_row],
      U238     = raw$U238[start_row:end_row],
      stringsAsFactors = FALSE
    )
    segment.data[, 2:5] <- lapply(segment.data[, 2:5],
                                  function(v) suppressWarnings(as.numeric(v)))
    ag <- counts_to_ages(segment.data$Pb206, segment.data$Pb207,
                         segment.data$U238, u238u235)
    segment.data$Age68 <- ag$Age68
    segment.data$Age75 <- ag$Age75
    segment.data$Age76 <- ag$Age76

  } else if ("Pb206_U238" %in% colnames(raw) &&
             !all(is.na(raw$Pb206_U238[start_row:end_row]))) {
    segment.data <- data.frame(
      Analysis   = raw$Analysis[start_row:end_row],
      Time       = raw$Time[start_row:end_row],
      Pb206U238  = raw$Pb206_U238[start_row:end_row],
      Pb207U235  = raw$Pb207_U235[start_row:end_row],
      Pb207Pb206 = raw$Pb207_Pb206[start_row:end_row],
      stringsAsFactors = FALSE
    )
    segment.data[, 2:5] <- lapply(segment.data[, 2:5],
                                  function(v) suppressWarnings(as.numeric(v)))
    ag <- ratios_to_ages(segment.data$Pb206U238, segment.data$Pb207U235,
                         segment.data$Pb207Pb206)
    segment.data$Age68 <- ag$Age68
    segment.data$Age75 <- ag$Age75
    segment.data$Age76 <- ag$Age76

  } else {
    stop("Unrecognized data format. Expected Age68, Pb206, or Pb206_U238 ",
         "columns; found: ", paste(colnames(raw), collapse = ", "),
         call. = FALSE)
  }

  extra_names <- setdiff(colnames(raw), colnames(segment.data))
  extra_data <- NULL
  if (length(extra_names) > 0) {
    extra_data <- raw[start_row:end_row, extra_names, drop = FALSE]
    extra_is_num <- vapply(extra_names, function(cn) {
      col_vals <- extra_data[[cn]]
      if (is.numeric(col_vals)) return(TRUE)
      converted <- suppressWarnings(as.numeric(as.character(col_vals)))
      if (length(converted) == 0) return(FALSE)
      sum(!is.na(converted)) / length(converted) > 0.5
    }, logical(1), USE.NAMES = FALSE)
    extra_names <- extra_names[extra_is_num]
  }

  list(data = segment.data, extra_names = extra_names, extra_data = extra_data)
}

# ---------------------------------------------------------------------------
#  Outlier detection
# ---------------------------------------------------------------------------

#' KPSS test statistic for level stationarity (with Bartlett window)
#' @keywords internal
kpss_stat <- function(y) {
  n <- length(y)
  e <- y - mean(y)
  S <- cumsum(e)
  l <- max(1L, floor(4 * (n / 100)^0.25))
  if (l >= n) l <- n - 1L
  if (l < 1L) l <- 1L
  g0 <- sum(e * e) / n
  s2 <- g0
  if (l >= 1L) {
    for (k in seq_len(l)) {
      gk <- sum(e[(k + 1L):n] * e[1:(n - k)]) / n
      s2 <- s2 + 2 * (1 - k / (l + 1)) * gk
    }
  }
  if (!is.finite(s2) || s2 <= 0) return(0)
  sum(S * S) / (n * n * s2)
}

#' Lightweight order selection for an ARIMA model
#'
#' Replicates the *intent* of `forecast::auto.arima()` (differencing order
#' from a KPSS test, then AIC search over p and q) without the dependency.
#' Grid search is used instead of the stepwise heuristic, so the selected
#' model can occasionally differ from auto.arima; for the purpose of residual
#' based outlier detection the effect is negligible.
#'
#' @param x Numeric vector.
#' @param max_p,max_q Maximum AR / MA order.
#' @param max_d Maximum differencing order.
#' @return Integer vector c(p, d, q).
#' @keywords internal
arima_order <- function(x, max_p = 2L, max_q = 2L, max_d = 1L) {
  d <- 0L
  y <- x
  repeat {
    if (length(y) < 12L || d >= max_d) break
    if (kpss_stat(y) < 0.463) break     # 5% critical value
    y <- diff(y)
    d <- d + 1L
  }

  best <- c(0L, d, 0L)
  best_aic <- Inf
  for (p in 0:max_p) {
    for (q in 0:max_q) {
      if (p == 0L && q == 0L) next
      fit <- suppressWarnings(try(stats::arima(x, order = c(p, d, q)),
                                  silent = TRUE))
      if (inherits(fit, "try-error")) next
      a <- fit$aic
      if (is.finite(a) && a < best_aic) {
        best_aic <- a
        best <- c(p, d, q)
      }
    }
  }
  best
}

#' ARIMA-based outlier detection
#'
#' Fits an automatically ordered ARIMA model, computes residual standard
#' deviations, and marks values beyond 2 SD as outliers (NA).
#'
#' @param series Numeric vector.
#' @param method \code{"arima"} (default, reproduces the published method),
#'   \code{"mad"} (fast robust z-score) or \code{"none"}.
#' @param n_sd Number of residual standard deviations used as the threshold.
#' @return Numeric vector with outliers replaced by NA.
#' @keywords internal
arima_outlier <- function(series, method = c("arima", "mad", "none"),
                          n_sd = 2) {
  method <- match.arg(method)
  if (identical(method, "none")) return(series)

  ok <- !is.null(series) && length(series) >= 10 &&
        sum(is.finite(series)) >= 10 &&
        is.finite(stats::sd(series, na.rm = TRUE)) &&
        stats::sd(series, na.rm = TRUE) > 0
  if (!ok) return(series)

  if (identical(method, "mad")) {
    med <- stats::median(series, na.rm = TRUE)
    mad <- stats::mad(series, na.rm = TRUE)
    if (!is.finite(mad) || mad == 0) return(series)
    series[abs(series - med) > n_sd * mad] <- NA
    return(series)
  }

  ord <- tryCatch(arima_order(series), error = function(e) c(0L, 0L, 0L))
  fit <- suppressWarnings(try(stats::arima(series, order = ord), silent = TRUE))
  if (inherits(fit, "try-error")) return(series)

  residuals <- suppressWarnings(try(stats::residuals(fit), silent = TRUE))
  if (inherits(residuals, "try-error") || length(residuals) != length(series)) {
    return(series)
  }

  res_sd <- stats::sd(residuals, na.rm = TRUE)
  if (!is.finite(res_sd) || res_sd == 0) return(series)

  # Align the residual vector with the original observations so that a fitted
  # ARIMA(p, d, q) still marks the right rows.
  len_res <- length(residuals)
  if (len_res < length(series)) {
    residuals <- c(rep(0, length(series) - len_res), residuals)
  }

  series[abs(residuals) > n_sd * res_sd] <- NA
  series
}

# ---------------------------------------------------------------------------
#  Preprocessing
# ---------------------------------------------------------------------------

#' Sanity-check a parsed segment table before it enters the pipeline
#'
#' Column-name typos, wrong time units and empty segments used to fail deep
#' inside the model fitting with unhelpful messages. This reports them up front
#' with something actionable.
#'
#' @param d Parsed data.frame with Time plus at least one age column.
#' @param label Sheet / chunk label used in the message.
#' @param strict If TRUE, stop on a hard error; if FALSE, return the problems
#'   as a character vector and let the caller warn and skip.
#' @param check_values Also test the age values themselves. Only meaningful
#'   once the ablation window has been applied - deep-profile files routinely
#'   contain nonsense ages in the gas-blank and wash-out portions.
#' @return Character vector of problems (empty when the data look sane).
#' @keywords internal
validate_segment_data <- function(d, label = "", strict = FALSE,
                                  check_values = FALSE) {
  p <- character(0)
  add <- function(...) p <<- c(p, paste0(...))

  n <- nrow(d)
  if (n < 10L) add("only ", n, " rows after parsing (need >= 10)")

  if (!"Time" %in% names(d)) {
    add("no Time column")
  } else {
    tt <- suppressWarnings(as.numeric(d$Time))
    if (all(is.na(tt))) {
      add("Time column is not numeric")
    } else {
      fin <- tt[is.finite(tt)]
      if (length(fin) >= 2L) {
        if (any(diff(fin) < 0)) {
          add("Time is not monotonically increasing - it may be in ",
              "milliseconds rather than seconds, or the rows are out of order")
        }
        rng <- range(fin)
        if (rng[2] > 1000) {
          add("Time spans ", signif(rng[2], 4),
              " - if that is milliseconds, divide by 1000 first")
        }
      }
    }
  }

  age_cols <- intersect(c("Age68", "Age75", "Age76", "Pb206", "Pb207",
                          "U238", "Pb206U238", "Pb207U235", "Pb207Pb206"),
                        names(d))
  if (length(age_cols) == 0L) {
    add("no recognisable age or isotope column")
  } else {
    for (cl in age_cols) {
      v <- suppressWarnings(as.numeric(d[[cl]]))
      if (all(is.na(v))) add(cl, " is entirely non-numeric")
    }
    ok <- rep(TRUE, n)
    for (cl in age_cols) ok <- ok & !is.na(suppressWarnings(as.numeric(d[[cl]])))
    if (sum(ok) < 10L) add("only ", sum(ok), " rows with a usable age value")

    if (isTRUE(check_values)) {
      for (cl in intersect(age_cols, c("Age68", "Age75", "Age76"))) {
        fin <- suppressWarnings(as.numeric(d[[cl]]))
        fin <- fin[is.finite(fin)]
        if (length(fin) == 0L) next
        neg <- sum(fin < 0)
        old <- sum(fin > 4568)
        if (neg > 0L) {
          add(neg, " negative ", cl,
              " value(s) inside the ablation window - check the common-Pb ",
              "correction")
        }
        if (old > 0L) {
          add(old, " ", cl, " value(s) above the age of the Earth inside the ",
              "ablation window - check the ratio columns")
        }
      }
    }
  }

  if (length(p) > 0L && isTRUE(strict)) {
    stop("Input problem in ", label, ": ", paste(p, collapse = "; "),
         call. = FALSE)
  }
  p
}

#' Discordance filter for U-Pb ages
#'
#' For ages < 1000 Ma: marks Age68 as NA if |Age68 - Age75| / Age75 > 0.1
#'
#' @param df data.frame with Age68, Age75, Age76 columns
#' @return Modified data.frame with `Age`, `subset_Age68`, `subset_Age76`
#' @keywords internal
discordance_filter <- function(df) {
  a75 <- df$Age75
  df$subset_Age68 <- ifelse(abs(df$Age68 - a75) / a75 <= 0.1, df$Age68, NA)
  df$subset_Age76 <- df$Age76
  df
}

#' Sliding window mean fill for NA values
#'
#' @param na_series Series with NA values to fill
#' @param original_series Original complete series for reference
#' @param window_size Half-window size (default 5)
#' @return Series with NAs filled by local mean
#' @keywords internal
mean_fill <- function(na_series, original_series, window_size = 5) {
  na_idx <- which(is.na(na_series))
  if (length(na_idx) == 0) return(na_series)
  local <- slide_mean_na(original_series, window_size)
  na_series[na_idx] <- local[na_idx]
  na_series
}

#' LOESS smoothing and standardisation
#'
#' @param df data.frame with subset_Age and Time columns
#' @param span LOESS span parameter (default 0.15)
#' @param family "gaussian" (ordinary least squares, the published pipeline) or
#'   "symmetric" (robust M-estimation, which downweights outliers). Using
#'   "symmetric" lets the ARIMA outlier screen be skipped entirely.
#' @return data.frame with loess_Age, standardized_loess, standardized_Age
#' @keywords internal
loess_segment <- function(df, span = 0.15, family = c("gaussian", "symmetric")) {
  family <- match.arg(family)
  loess_model <- try(
    stats::loess(subset_Age ~ Time, data = df, span = span, family = family),
    silent = TRUE)
  if (inherits(loess_model, "try-error")) {
    warning("LOESS smoothing failed; falling back to raw ages.", call. = FALSE)
    df$loess_Age <- df$subset_Age
  } else {
    df$loess_Age <- stats::predict(loess_model)
  }
  df <- df[complete.cases(df$loess_Age), ]
  if (nrow(df) == 0) return(df)

  df$log_loess <- df$loess_Age
  df$log_Age   <- df$subset_Age

  minloess <- min(df$log_loess)
  maxloess <- max(df$log_loess)
  rng <- maxloess - minloess

  if (!is.finite(rng) || rng <= 0) {
    df$standardized_loess <- 0
    df$standardized_Age   <- 0
  } else {
    df$standardized_loess <- (df$log_loess - minloess) / rng
    df$standardized_Age   <- (df$log_Age   - minloess) / rng
  }
  df <- df[complete.cases(df$standardized_loess), ]

  attr(df, "minloess") <- minloess
  attr(df, "maxloess") <- maxloess
  df
}

# ---------------------------------------------------------------------------
#  Segmentation (vectorised)
# ---------------------------------------------------------------------------

#' Build plateau segment definitions from changepoints
#'
#' @param df data.frame with Time and standardized_loess columns
#' @param changepoints Integer vector of changepoint row indices
#' @return list(segments, df, segment_starts, segment_ends)
#' @keywords internal
build_segments <- function(df, changepoints) {
  segment_starts <- c(1L, utils::head(changepoints, -1L) + 1L)
  segment_ends   <- changepoints
  m <- length(segment_starts)

  segments <- data.frame(
    Start = df$Time[segment_starts],
    End   = df$Time[segment_ends],
    stringsAsFactors = FALSE
  )

  segments$standardized_Mean <- seg_mean(df$standardized_loess,
                                         segment_starts, segment_ends)

  df$IS_loess <- seg_normalise(df$loess_Age, segment_starts, segment_ends)

  minloess <- attr(df, "minloess")
  maxloess <- attr(df, "maxloess")
  df$Restored_loess <- df$standardized_loess * (maxloess - minloess) + minloess
  df$Restored_age   <- df$standardized_Age   * (maxloess - minloess) + minloess

  segments$Segment_Mean <- seg_mean(df$Restored_loess,
                                    segment_starts, segment_ends)

  segments$Number <- seq_len(m)
  segments$Time_step <- segments$End - segments$Start
  segments$Max_step  <- max(segments$Time_step)
  segments$Min_step  <- min(segments$Time_step)

  list(segments = segments, df = df, segment_starts = segment_starts,
       segment_ends = segment_ends)
}

#' Slope and intercept for each plateau segment
#'
#' Closed-form OLS per segment, vectorised (no `lm()` call per segment).
#'
#' @keywords internal
calc_slopes <- function(df, segments, seg_starts = NULL, seg_ends = NULL) {
  if (is.null(seg_starts)) {
    seg_starts <- match(segments$Start, df$Time)
    seg_ends   <- match(segments$End, df$Time)
  }
  a <- seg_lm(df$Time, df$standardized_loess, seg_starts, seg_ends)
  segments$standardized_intercept <- a$intercept
  segments$standardized_slope     <- a$slope

  b <- seg_lm(df$Time, df$loess_Age, seg_starts, seg_ends)
  segments$loess_intercept <- b$intercept
  segments$loess_slope     <- b$slope
  segments
}

#' Plateau variance (within-segment normalised LOESS)
#' @keywords internal
calc_variance <- function(df, seg_starts, seg_ends) {
  seg_var(df$IS_loess, seg_starts, seg_ends)
}

#' Plateau uncertainties
#'
#' The three columns that existed before are kept with exactly the same
#' meaning and values:
#'   Calibration uncertainty = Segment_Mean * calibration_uncertainty
#'   Plateau uncertainty     = sd(Raw_Age) / sqrt(n) inside the plateau
#'   Total uncertainty       = sqrt(cal^2 + plateau^2)
#'
#' v1.2.0 adds a decomposition separating what was *measured* from what was
#' *assumed*, plus the systematic contribution of the decay constants:
#'
#'   Random uncertainty      = sqrt(plateau^2 + cal^2)
#'                             (use this to compare ages within a session)
#'   Systematic uncertainty  = age * relative 1-sigma of the decay constants
#'                             (fully correlated between samples)
#'   Total incl. decay       = sqrt(random^2 + systematic^2)
#'
#' The legacy `Total_uncertainty` deliberately excludes the decay constants so
#' that previously published numbers do not move. Quote the new column
#' (`Total_uncertainty_full`) in a paper.
#'
#' @param segments Segment data.frame.
#' @param df Point-level data.frame (needs Raw_Age, Age68, Age76).
#' @param seg_starts,seg_ends Integer segment boundaries.
#' @param calibration_uncertainty Relative 1-sigma reproducibility of the
#'   primary reference material (default 0.03, i.e. 3 \%).
#' @param decay_system Age system behind the plateau age: "auto", "Age68",
#'   "Age75" or "Age76". "auto" picks Age68 below 1 Ga and Age76 above it,
#'   matching how Raw_Age is built.
#' @keywords internal
calc_uncertainty <- function(segments, df, seg_starts, seg_ends,
                             calibration_uncertainty = 0.03,
                             decay_system = c("auto", "Age68", "Age75",
                                              "Age76")) {
  decay_system <- match.arg(decay_system)

  segments$Calibration_uncertainty <-
    segments$Segment_Mean * calibration_uncertainty
  n  <- seg_count(df$Raw_Age, seg_starts, seg_ends)
  sd <- seg_sd_na(df$Raw_Age, seg_starts, seg_ends)
  pu <- sd / sqrt(n)
  pu[n < 2L] <- NA_real_
  segments$Plateau_uncertainty <- pu
  segments$Total_uncertainty <-
    sqrt(segments$Calibration_uncertainty^2 + pu^2)

  # ---- v1.2.0 decomposition ------------------------------------------------
  m <- nrow(segments)
  if (identical(decay_system, "auto")) {
    seg_id_vec <- rep(NA_integer_, nrow(df))
    for (i in seq_along(seg_starts)) {
      seg_id_vec[seg_starts[i]:seg_ends[i]] <- i
    }
    is68 <- !is.na(df$Age68) & df$Age68 < 1000
    frac68 <- tapply(is68, seg_id_vec, mean, na.rm = TRUE)
    frac68 <- as.numeric(frac68[as.character(seq_len(m))])
    system <- ifelse(is.na(frac68) | frac68 >= 0.5, "Age68", "Age76")
  } else {
    system <- rep(decay_system, m)
  }

  rel <- adept_decay_rel_1s(system, segments$Segment_Mean)
  segments$Decay_system <- as.character(system)
  segments$Decay_constant_uncertainty <- segments$Segment_Mean * rel
  segments$Random_uncertainty <- segments$Total_uncertainty
  segments$Systematic_uncertainty <- segments$Decay_constant_uncertainty
  segments$Total_uncertainty_full <-
    sqrt(segments$Random_uncertainty^2 + segments$Systematic_uncertainty^2)
  segments$Relative_uncertainty_pct <-
    100 * segments$Total_uncertainty_full / segments$Segment_Mean
  segments
}

#' Means of extra (non-age) columns per plateau segment
#' @keywords internal
calc_extra_means <- function(segments, df, extra_names,
                             seg_starts = NULL, seg_ends = NULL) {
  if (length(extra_names) == 0) return(segments)
  if (is.null(seg_starts)) {
    seg_starts <- match(segments$Start, df$Time)
    seg_ends   <- match(segments$End, df$Time)
  }
  for (col_name in extra_names) {
    mean_name <- paste0(col_name, "_Mean")
    if (col_name %in% colnames(df)) {
      segments[[mean_name]] <-
        seg_mean_na(df[[col_name]], seg_starts, seg_ends)$mean
    } else {
      segments[[mean_name]] <- NA_real_
    }
  }
  segments
}

#' Mean and total uncertainty of Age68/Age75/Age76 per plateau
#' @keywords internal
calc_age_means <- function(segments, df, seg_starts = NULL, seg_ends = NULL) {
  if (is.null(seg_starts)) {
    seg_starts <- match(segments$Start, df$Time)
    seg_ends   <- match(segments$End, df$Time)
  }
  for (age_col in c("Age68", "Age75", "Age76")) {
    if (!age_col %in% colnames(df)) next
    v  <- df[[age_col]]
    mn <- seg_mean_na(v, seg_starts, seg_ends)
    n  <- mn$n
    s  <- seg_sd_na(v, seg_starts, seg_ends)
    plateau_unc <- s / sqrt(n)
    cal_unc     <- mn$mean * 0.03
    unc <- sqrt(cal_unc^2 + plateau_unc^2)
    unc[n < 2L] <- NA_real_
    segments[[paste0(age_col, "_Mean")]] <- mn$mean
    segments[[paste0(age_col, "_Total_uncertainty")]] <- unc
  }
  segments
}

#' U-Pb concordance (Age68 / Age75 * 100) per plateau
#' @keywords internal
calc_concordance <- function(segments) {
  segments$Concordance <- (segments$Age68_Mean / segments$Age75_Mean) * 100
  segments
}

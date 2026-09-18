# Tests for the per-zircon pipeline steps (zircon.R).
#
# `adept()` was a single 700-line function whose body interleaved bookkeeping
# with four different kinds of work. These tests target the seams that were
# cut out of it, which the end-to-end tests in test-adept.R cannot reach: the
# exact reason a zircon is dropped, and the grouping rule that decides what a
# "zircon" is in the first place.
#
# The configuration list is built here rather than by calling adept(), so that
# a step can be exercised without a file on disk.

cfg_for <- function(...) {
  base <- list(
    u238u235                = ADEPT_U238U235,
    validate_input          = FALSE,
    lower_ablation_time     = 0,
    upper_ablation_time     = Inf,
    preprocess              = "arima_loess",
    smooth                  = "loess",
    outlier_method          = "arima",
    outlier_sd              = 2,
    calibration_uncertainty = 0.03,
    min_age_limit           = 0,
    max_age_limit           = 4540,
    variance_threshold      = 0.1192,
    min_plateau_resolution  = NULL,
    filter_direction        = "Forward",
    direction_method        = "monotonic",
    direction_tolerance     = 0.02,
    all_extra_names         = character(0)
  )
  utils::modifyList(base, list(...))
}

zircon_frame <- function(age68, time = seq_along(age68)) {
  d <- data.frame(Analysis = rep("A", length(age68)), Time = time,
                  Age68 = age68, stringsAsFactors = FALSE)
  d$.ROWID. <- seq_len(nrow(d))
  d
}

# ---------------------------------------------------------------------------
#  Grouping
# ---------------------------------------------------------------------------

test_that("a sheet is split by Analysis, never by row count", {
  # 11 rows, three zircons, chunk_size smaller than two of them. Cutting every
  # chunk_size rows used to merge different zircons and split one in half.
  d <- data.frame(Analysis = rep(c("A", "B", "C"), c(4, 2, 5)), Time = 1:11)
  g <- adept_zircon_groups(d, chunk_size = 3L, "Sheet1")

  expect_equal(g$n, 3L)
  expect_equal(g$starts, c(1L, 5L, 7L))
  expect_equal(g$ends,   c(4L, 6L, 11L))
  expect_equal(g$labels, c("A", "B", "C"))
  expect_equal(g$analysis_col, "Analysis")
})

test_that("consecutive rows repeat an Analysis value; a later return is a new zircon", {
  # rle() semantics, made explicit: the same sample name appearing again after
  # another one is two zircons, not one.
  d <- data.frame(Analysis = c("A", "A", "B", "A"), Time = 1:4)
  g <- adept_zircon_groups(d, 411L, "S")
  expect_equal(g$n, 3L)
  expect_equal(g$labels, c("A", "B", "A"))
})

test_that("chunk_size is the fallback only when there is no Analysis column", {
  d <- data.frame(Time = 1:10)
  g <- adept_zircon_groups(d, chunk_size = 4L, "S1")

  expect_true(is.na(g$analysis_col))
  expect_equal(g$n, 3L)
  expect_equal(g$starts, c(1L, 5L, 9L))
  expect_equal(g$ends,   c(4L, 8L, 10L))
  expect_equal(g$labels, c("S1_1", "S1_5", "S1_9"))
})

test_that("'Analysis_' is accepted as an alias for 'Analysis'", {
  d <- data.frame(Analysis_ = c("A", "A", "B"), Time = 1:3)
  g <- adept_zircon_groups(d, 411L, "S")
  expect_equal(g$n, 2L)
  expect_equal(g$analysis_col, "Analysis_")
})

test_that("a degenerate sheet still yields one coherent group", {
  # An empty sheet has no runs at all. `starts` and `ends` must stay the same
  # length, or the main loop would index past the end of one of them.
  d <- data.frame(Analysis = character(0), Time = numeric(0))
  g <- adept_zircon_groups(d, 411L, "S")
  expect_equal(g$n, 1L)
  expect_equal(length(g$starts), length(g$ends))
  expect_equal(g$ends, 0L)

  # The row-count fallback has to agree: neither branch may report zero units,
  # or the progress denominator and the loop disagree.
  e <- adept_zircon_groups(data.frame(Time = numeric(0)), 411L, "S")
  expect_equal(e$n, 1L)
  expect_equal(e$ends, 0L)
})

# ---------------------------------------------------------------------------
#  The "no ages" guard
# ---------------------------------------------------------------------------

test_that("absent age columns are not tested for NA", {
  # all(is.na(NULL)) is TRUE. Testing an absent column therefore rejected a
  # perfectly good input, which is what made a two-column Format 4 sheet
  # unreadable.
  d <- data.frame(Age68 = c(1, 2), Age75 = c(NA, NA))
  expect_true(adept_has_ages(d, has_75 = FALSE, has_76 = FALSE))
  expect_false(adept_has_ages(d, has_75 = TRUE, has_76 = FALSE))

  # A present column with a value counts, even if the column is mostly blank.
  expect_true(adept_has_ages(data.frame(Age68 = c(1, NA), Age75 = c(2, NA)),
                             has_75 = TRUE, has_76 = FALSE))

  expect_false(adept_has_ages(data.frame(Age68 = c(NA_real_, NA_real_)),
                              has_75 = FALSE, has_76 = FALSE))
})

# ---------------------------------------------------------------------------
#  Ablation window and Raw_Age
# ---------------------------------------------------------------------------

test_that("the window trims the profile and Raw_Age uses Age68 below 1 Ga", {
  d <- zircon_frame(rep(100, 6), time = c(1, 5, 29, 40, 58, 70))
  d$Age75 <- 100; d$Age76 <- 100
  p <- adept_prepare_window(d, cfg_for(lower_ablation_time = 29,
                                       upper_ablation_time = 58),
                            "S", 1L, 6L)

  expect_true(p$ok)
  expect_equal(nrow(p$data), 3L)                    # 29, 40, 58
  expect_equal(p$data$Time, c(29, 40, 58))
  expect_equal(p$data$Raw_Age, rep(100, 3))
  expect_true(p$has_75); expect_true(p$has_76)
})

test_that("Raw_Age picks Age76 above 1 Ga and mirrors the branch for its sigma", {
  # A Mesoarchaean grain: 206Pb/238U is too imprecise, so the age comes from
  # 207Pb/206Pb. The sigma has to follow whichever branch was taken.
  d <- zircon_frame(c(100, 2000, 2000, 100, 100))
  d$Age68_1s <- c(1, 2, 3, 4, 5)
  d$Age76    <- c(2000, 2000, 100, 2000, 2000)
  d$Age76_1s <- c(10, 20, 30, 40, 50)

  p <- adept_prepare_window(d, cfg_for(), "S", 1L, 5L)
  expect_true(p$ok)
  expect_equal(p$data$Raw_Age, c(100, 2000, NA, 100, 100))
  expect_equal(p$data$Raw_Age_sigma, c(1, 20, NA, 4, 5))
})

test_that("a pre-reduced input with no 207Pb ages enters the pipeline", {
  d <- zircon_frame(rep(300, 5))
  d$Age68_1s <- rep(2, 5)

  p <- adept_prepare_window(d, cfg_for(), "S", 1L, 5L)
  expect_true(p$ok)
  expect_false(p$has_75)
  expect_false(p$has_76)
  expect_equal(p$data$Raw_Age_sigma, rep(2, 5))
})

test_that("a window that removes every row says so", {
  d <- zircon_frame(rep(100, 5), time = 1:5)
  d$Age75 <- 100; d$Age76 <- 100
  p <- adept_prepare_window(d, cfg_for(lower_ablation_time = 100,
                                       upper_ablation_time = 200),
                            "S", 1L, 5L)
  expect_false(p$ok)
  expect_equal(p$label, "empty window")
  expect_equal(p$points, 0L)
})

test_that("an all-blank Age68 column is reported before any fitting", {
  d <- zircon_frame(rep(NA_real_, 20))
  p <- adept_prepare_window(d, cfg_for(), "S", 1L, 20L)
  expect_false(p$ok)
  expect_equal(p$label, "no ages")
})

# ---------------------------------------------------------------------------
#  Whole-zircon entry point
# ---------------------------------------------------------------------------

test_that("one zircon in, one segment table out", {
  # Two flat domains, 60 points each, with the sigma that makes MSWD defined.
  raw <- data.frame(Analysis = rep("A", 120), Time = seq(0, 59.5, length.out = 120))
  raw$Age68 <- c(rep(100, 60), rep(140, 60))
  raw$Age75 <- raw$Age68
  raw$Age76 <- raw$Age68
  raw$Age68_1s <- 2; raw$Age75_1s <- 2; raw$Age76_1s <- 2

  r <- adept_one_zircon(raw, 1L, 120L, "S", cfg_for(smooth = "none"))
  expect_true(r$ok)
  expect_equal(r$analysis, "A")
  expect_equal(nrow(r$data), 120L)
  expect_true(all(c("MSWD", "MSWD_prob", "Uncertainty_n") %in%
                    names(r$segments)))

  # PELT splits this into three: the two domains, plus the single point sitting
  # on the discontinuity, which the ARIMA screen blanked and mean_fill() then
  # put back. The one-point segment is dropped by the duration and variance
  # steps (Time_step 0, variance undefined for a single point), so two
  # plateaus survive — and their ages are the two domain ages exactly.
  confirmed <- r$segments[!is.na(r$segments$Filter_4), ]
  expect_equal(nrow(confirmed), 2L)
  expect_equal(confirmed$Segment_Mean, c(100, 140), tolerance = 1e-9)
})

test_that("a perfectly flat segment has variance zero, not a missing value", {
  # seg_normalise() returns NaN for a zero-range segment, which propagates
  # through `ifelse(Variance <= threshold, ...)` as NA and silently discards
  # the segment. A flat domain is the best possible plateau, so it must score
  # zero. Nothing on a LOESS-smoothed profile is ever exactly flat; this is the
  # shape a down-hole-fractionation corrected domain has under smooth = "none".
  #
  # Reproduces the two lines build_segments() runs, in order.
  x      <- c(rep(5, 30), rep(9, 30), 7)
  starts <- c(1L, 31L, 61L)
  ends   <- c(30L, 60L, 61L)
  df     <- data.frame(loess_Age = x)
  df$IS_loess <- seg_normalise(df$loess_Age, starts, ends)

  expect_true(all(is.nan(df$IS_loess[1:60])))            # the cause
  expect_true(is.na(seg_var(df$IS_loess, starts, ends)[1]))

  v <- calc_variance(df, starts, ends)
  expect_equal(v[1:2], c(0, 0))                          # the fix
  # A single-point segment is a different case: variance is genuinely
  # undefined there and must stay missing rather than become 0.
  expect_true(is.na(v[3]))
})

test_that("a sheet in no recognised format is reported, not thrown", {
  raw <- data.frame(Analysis = "A", Time = 1:20, Foo = 1:20)
  expect_warning(r <- adept_one_zircon(raw, 1L, 20L, "S", cfg_for()),
                 "skipped")
  expect_false(r$ok)
  expect_equal(r$label, "skipped")
  expect_equal(r$points, 0L)
})

test_that("fewer than ten usable ages is reported as such", {
  # Only five rows carry an age. The guard that catches this had been written
  # twice in adept(), and fixing one copy left the other in place.
  raw <- data.frame(Analysis = rep("A", 20), Time = 1:20,
                    Age68 = c(rep(100, 5), rep(NA_real_, 15)),
                    Age75 = 100, Age76 = 100)
  r <- adept_one_zircon(raw, 1L, 20L, "S", cfg_for(smooth = "none"))
  expect_false(r$ok)
  expect_equal(r$label, "too few points")
  expect_equal(r$points, 5L)
})

test_that("every failure label is a short, stable string", {
  # These become the progress messages, so they are part of the interface.
  raw <- data.frame(Analysis = rep("A", 12), Time = 1:12, Age68 = rep(100, 12),
                    Age75 = rep(100, 12), Age76 = rep(100, 12))
  r <- adept_one_zircon(raw, 1L, 12L, "S",
                        cfg_for(lower_ablation_time = 99,
                                upper_ablation_time = 100))
  expect_false(r$ok)
  expect_equal(r$label, "empty window")
  expect_true(is.character(r$label) && nchar(r$label) < 32)
})

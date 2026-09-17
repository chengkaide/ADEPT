# End-to-end regression tests on the bundled example workbook.
#
# The headline numbers below are the ones the package has produced since
# v1.1.0. They are the anchor that stops a refactor from silently changing
# somebody's published age.

example_file <- function() {
  f <- system.file("extdata", "Input(1sample).xlsx", package = "ADEPT")
  if (f == "") skip("example workbook not installed")
  f
}

run_example <- function(make_plots = FALSE, ...) {
  suppressWarnings(adept(example_file(), output_path = NA,
                         make_plots = make_plots, verbose = FALSE, ...))
}

test_that("the default run reproduces the published ages", {
  r <- run_example()
  expect_equal(nrow(r$summary), 3)
  expect_equal(r$summary[["Final age (Ma)"]],
               c(21.61724, 22.56383, 23.60464), tolerance = 1e-5)
  expect_equal(r$summary[["Final total uncertainty (Ma)"]],
               c(0.6886189, 0.6980347, 0.7447352), tolerance = 1e-6)
  expect_equal(r$summary$Analysis, rep("DJ27_18", 3))
  expect_equal(r$summary$Points, rep(161, 3))
})

test_that("numeric output columns really are numeric", {
  # An earlier version built each row with c() and turned every number into
  # a string, which made the exported workbook useless for further work.
  r <- run_example()
  num_cols <- setdiff(colnames(r$full), c("Analysis", "Decay system"))
  for (cn in num_cols) {
    expect_true(is.numeric(r$full[[cn]]), info = cn)
  }
  expect_true(is.character(r$full[["Decay system"]]))
})

test_that("the uncertainty budget decomposes as documented", {
  r <- run_example()
  full <- r$full
  keep <- !is.na(full[["Final age (Ma)"]])

  random <- full[["Random uncertainty (Ma)"]][keep]
  system <- full[["Systematic uncertainty (Ma)"]][keep]
  total  <- full[["Final total uncertainty incl. decay (Ma)"]][keep]
  legacy <- full[["Final total uncertainty (Ma)"]][keep]
  decay  <- full[["Decay constant uncertainty (Ma)"]][keep]

  expect_equal(random, legacy)                       # random == the old total
  expect_equal(system, decay)                        # only decay is systematic
  expect_equal(total, sqrt(random^2 + system^2))
  expect_true(all(total > legacy))                   # adding decay only grows it
  expect_true(all(decay > 0))
  expect_equal(unique(full$`Decay system`[keep]), "Age68")
  expect_equal(full[["Relative uncertainty (%)"]][keep],
               100 * total / full[["Final age (Ma)"]][keep])
})

test_that("calibration_uncertainty is honoured and defaults to 3 percent", {
  r3 <- run_example()
  r1 <- run_example(calibration_uncertainty = 0.015)

  age3 <- r3$summary[["Final age (Ma)"]]
  age1 <- r1$summary[["Final age (Ma)"]]
  expect_equal(age3, age1)                            # ages do not depend on it

  # halving the assumed reproducibility nearly halves the random term
  rat <- r1$summary[["Random uncertainty (Ma)"]] /
         r3$summary[["Random uncertainty (Ma)"]]
  expect_true(all(rat > 0.5 & rat < 0.75))

  expect_error(run_example(calibration_uncertainty = -1), "non-negative")
  expect_error(run_example(calibration_uncertainty = c(1, 2)), "single")
})

test_that("direction_method = strict reproduces the old selector", {
  a <- run_example(direction_method = "monotonic")
  b <- run_example(direction_method = "strict")
  expect_equal(a$summary[["Final age (Ma)"]],
               b$summary[["Final age (Ma)"]])
})

test_that("Reverse returns a subset of the Forward plateaus here", {
  fwd <- run_example(filter_direction = "Forward")
  rev_ <- run_example(filter_direction = "Reverse")
  expect_true(nrow(rev_$summary) >= 1)
  expect_true(nrow(rev_$summary) <= nrow(fwd$summary))
})

test_that("robust_loess runs and stays geologically plausible", {
  r <- run_example(preprocess = "robust_loess")
  expect_true(nrow(r$summary) >= 1)
  ages <- r$summary[["Final age (Ma)"]]
  expect_true(all(ages > 15 & ages < 30),
              info = paste(ages, collapse = ", "))
})

test_that("the output workbook is written and read back", {
  f <- tempfile(fileext = ".xlsx")
  on.exit(unlink(f), add = TRUE)
  suppressWarnings(adept(example_file(), output_path = f,
                         make_plots = FALSE, verbose = FALSE))
  expect_true(file.exists(f))

  back <- adept_read(f)
  expect_equal(names(back), c("Summary", "Full_Results"))
  expect_equal(back$Summary[["Final age (Ma)"]],
               c(21.61724, 22.56383, 23.60464), tolerance = 1e-5)
  expect_true(is.numeric(back$Full_Results[["Final age (Ma)"]]))
  # every v1.2.0 column survives the round trip, including the long names
  expect_true(all(c("Final total uncertainty incl. decay (Ma)",
                    "Random uncertainty (Ma)", "Systematic uncertainty (Ma)",
                    "Decay constant uncertainty (Ma)", "Decay system",
                    "Relative uncertainty (%)", "Confirmed plateaus",
                    "Total segments") %in% colnames(back$Full_Results)))
})

test_that("only confirmed plateaus appear in the summary sheet", {
  r <- run_example()
  expect_true(all(!is.na(r$summary[["Final age (Ma)"]])))
  expect_equal(nrow(r$summary),
               sum(!is.na(r$full[["Final age (Ma)"]])))
})

test_that("plots are lightweight records that can be drawn", {
  r <- run_example(make_plots = TRUE)
  expect_equal(length(r$plots), 1L)
  expect_s3_class(r$plots[[1]], "adept_profile")
  f <- tempfile(fileext = ".pdf")
  on.exit(unlink(f), add = TRUE)
  grDevices::pdf(f)
  draw_profile(r$plots[[1]])
  grDevices::dev.off()
  expect_true(file.size(f) > 0)
})

test_that("validate_segment_data flags the problems it is meant to", {
  n <- 200
  bad <- data.frame(Analysis = "X",
                    Time = seq(0, 82000, length.out = n),
                    Age68 = rnorm(n, 800, 20),
                    Age75 = rnorm(n, 810, 20),
                    Age76 = rnorm(n, 9000, 20))
  p <- validate_segment_data(bad, "bad", check_values = TRUE)
  expect_true(any(grepl("milliseconds", p)))
  expect_true(any(grepl("age of the Earth", p)))

  # values are only checked when asked, because wash-out rows are nonsense
  p2 <- validate_segment_data(bad, "bad")
  expect_true(any(grepl("milliseconds", p2)))
  expect_false(any(grepl("age of the Earth", p2)))

  good <- data.frame(Analysis = "Y", Time = seq(0, 80, length.out = 100),
                     Age68 = rnorm(100, 800, 20),
                     Age75 = rnorm(100, 810, 20),
                     Age76 = rnorm(100, 820, 20))
  expect_length(validate_segment_data(good, "good", check_values = TRUE), 0)
  expect_error(validate_segment_data(good[0, ], "empty", strict = TRUE),
               "Input problem")
})

test_that("adept_sensitivity returns one row per combination", {
  s <- adept_sensitivity(example_file(),
                         variance_threshold = c(0.1192, 0.25),
                         min_plateau_resolution = 5,
                         filter_direction = "Forward",
                         verbose = FALSE)
  expect_equal(nrow(s), 2)
  expect_true(all(c("variance_threshold", "n_plateaus", "ages", "mean_age") %in%
                    colnames(s)))
  expect_true(all(s$n_plateaus > 0))
  # a looser threshold cannot yield fewer plateaus
  expect_true(s$n_plateaus[2] >= s$n_plateaus[1])
})

test_that("missing files and bad parameters fail loudly", {
  expect_error(adept(tempfile(fileext = ".xlsx"), output_path = NA,
                     make_plots = FALSE, verbose = FALSE), "not found")
  expect_error(run_example(u238u235 = -1), "positive")
})

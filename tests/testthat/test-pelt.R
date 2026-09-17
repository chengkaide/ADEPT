# Regression tests for the built-in PELT implementation.
#
# The expected changepoint vectors below were produced by this implementation
# and independently confirmed against changepoint::cpt.mean(method = "PELT",
# penalty = "Manual"). They are hard-coded so the test suite does not need
# the changepoint package installed.

step_series <- function() {
  set.seed(123)
  x <- rnorm(100)
  x[41:70] <- x[41:70] + 3
  x
}

three_level_series <- function() {
  set.seed(7)
  c(rnorm(60, 0, 1), rnorm(60, 4, 1), rnorm(60, -2, 1))
}

test_that("changepoints are reproduced for a known single step", {
  x <- step_series()
  pen <- 0.1 * log(100)

  expect_equal(pelt_mean(x, pen = pen, minseglen = 1L),
               c(2L, 3L, 5L, 6L, 7L, 10L, 11L, 14L, 15L, 16L, 17L, 18L, 19L,
                 25L, 26L, 28L, 29L, 30L, 32L, 37L, 40L, 43L, 44L, 45L, 48L,
                 53L, 54L, 55L, 56L, 57L, 61L, 65L, 69L, 70L, 71L, 72L, 73L,
                 75L, 76L, 78L, 86L, 88L, 89L, 93L, 94L, 95L, 96L, 98L,
                 100L))
  expect_equal(pelt_mean(x, pen = pen, minseglen = 5L),
               c(7L, 12L, 17L, 26L, 40L, 48L, 56L, 65L, 70L, 75L, 85L, 100L))
  expect_equal(pelt_mean(x, pen = pen, minseglen = 10L),
               c(17L, 29L, 40L, 56L, 70L, 83L, 100L))
})

test_that("changepoints are reproduced for a three-level series", {
  y <- three_level_series()
  pen <- 0.5 * log(180)

  expect_equal(pelt_mean(y, pen = pen, minseglen = 5L),
               c(9L, 15L, 60L, 71L, 78L, 120L, 180L))
  expect_equal(pelt_mean(y, pen = pen, minseglen = 20L),
               c(28L, 60L, 120L, 180L))
  expect_equal(pelt_mean(y, pen = pen, minseglen = 61L),
               c(119L, 180L))
})

test_that("every segment respects minseglen", {
  # This is the contract the v1.1.0 implementation broke: it could emit a
  # first segment shorter than minseglen.
  set.seed(99)
  for (trial in seq_len(60)) {
    n  <- sample(c(10L, 20L, 30L, 45L, 90L), 1)
    ms <- sample(seq_len(max(1L, n %/% 2L)), 1)
    x  <- rnorm(n)
    cp <- pelt_mean(x, pen = 0.04 * log(n), minseglen = ms)
    lens <- diff(c(0L, cp))
    expect_true(all(lens >= ms),
                info = sprintf("n=%d minseglen=%d lengths=%s", n, ms,
                               paste(lens, collapse = ",")))
    expect_equal(tail(cp, 1), n)
    expect_true(all(cp >= 1L & cp <= n))
  }
})

test_that("a zero or negative penalty is rejected", {
  # With pen = 0 the pruning step can discard the only feasible candidate,
  # so it is refused rather than silently returning a wrong answer.
  expect_error(pelt_mean(rnorm(20), pen = 0), "must be a single positive")
  expect_error(pelt_mean(rnorm(20), pen = -1), "must be a single positive")
  expect_error(pelt_mean(rnorm(20), pen = NA_real_), "must be a single positive")
})

test_that("the result is deterministic", {
  x <- step_series()
  expect_identical(pelt_mean(x, pen = 0.5),
                   pelt_mean(x, pen = 0.5))
})

test_that("edge cases do not error", {
  expect_equal(pelt_mean(5, pen = 0.1), 1L)          # single point
  # the return value is the vector of segment ends, not a count
  expect_equal(pelt_mean(c(1, 5), pen = 0.1), c(1L, 2L))
  # minseglen larger than the series: one segment
  expect_equal(pelt_mean(rnorm(10), pen = 0.1, minseglen = 20L), 10L)
})

test_that("pelt_segmentation uses the SAIC penalty and degrades gracefully", {
  x <- step_series()
  s <- pelt_segmentation(x, length(x))
  expect_equal(s$SAIC, 0.04 * log(length(x)))
  expect_equal(tail(s$changepoints, 1), length(x))

  # constant input has no changepoints to find
  expect_equal(pelt_segmentation(rep(1, 50), 50)$changepoints, 50L)
})

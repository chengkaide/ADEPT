# Tests for the vectorised segmented statistics.
#
# Every function here replaces an sapply/lm loop, so the tests compare against
# the naive implementation they were written to replace.

starts_ends <- function(n = 411, n_seg = 40) {
  starts <- round(seq(1, n, length.out = n_seg + 1L))
  list(starts = starts[-length(starts)],
       ends   = starts[-1] - 1L)
}

test_that("seg_mean, seg_var and seg_sd match the naive versions", {
  set.seed(1)
  x <- rnorm(411)
  se <- starts_ends()
  s <- se$starts; e <- se$ends

  expect_equal(seg_mean(x, s, e),
               vapply(seq_along(s), function(i) mean(x[s[i]:e[i]]), numeric(1)))
  expect_equal(seg_var(x, s, e),
               vapply(seq_along(s), function(i) var(x[s[i]:e[i]]), numeric(1)))
  expect_equal(seg_sd_na(x, s, e),
               vapply(seq_along(s), function(i) sd(x[s[i]:e[i]]), numeric(1)))
})

test_that("seg_lm matches lm() including intercept", {
  set.seed(2)
  x <- rnorm(411)
  tt <- seq(0, 82, length.out = 411)
  se <- starts_ends()
  s <- se$starts; e <- se$ends

  co <- vapply(seq_along(s), function(i)
    coef(lm(x[s[i]:e[i]] ~ tt[s[i]:e[i]])), numeric(2))
  got <- seg_lm(tt, x, s, e)

  expect_equal(got$slope, co[2, ], tolerance = 1e-6)
  expect_equal(got$intercept, co[1, ], tolerance = 1e-6)
})

test_that("single-point and flat segments give NA slopes like lm()", {
  tt <- seq(0, 82, length.out = 20)
  y  <- rnorm(20)
  y[10:15] <- 5                      # segment 2 is genuinely flat
  # segment 1 is a single point, segment 2 is flat
  got <- seg_lm(tt, y, c(5L, 10L), c(5L, 15L))
  expect_true(is.na(got$slope[1]))
  expect_true(is.na(got$intercept[1]))
  expect_equal(got$slope[2], 0, tolerance = 1e-12)
  expect_equal(got$intercept[2], 5, tolerance = 1e-12)
})

test_that("seg_normalise reproduces the lapply it replaces", {
  set.seed(3)
  x <- rnorm(60)
  s <- c(1L, 21L, 41L); e <- c(20L, 40L, 60L)

  old <- unlist(lapply(seq_along(s), function(i) {
    seg <- x[s[i]:e[i]]
    (seg - min(seg)) / (max(seg) - min(seg))
  }))
  expect_equal(seg_normalise(x, s, e), old)
  expect_true(all(seg_normalise(x, s, e) >= 0))
  expect_true(all(seg_normalise(x, s, e) <= 1))
})

test_that("slide_mean_na matches the sliding-window mean", {
  set.seed(4)
  x <- rnorm(50)
  x[c(5, 6, 20, 49)] <- NA
  k <- 5L

  got <- slide_mean_na(x, k)
  n <- length(x)
  want <- vapply(seq_len(n), function(i)
    mean(x[max(1, i - k):min(n, i + k)], na.rm = TRUE), numeric(1))

  expect_equal(got, want)
})

test_that("seg_count and seg_mean_na ignore NA", {
  x <- c(1, NA, 3, 4, NA, 6)
  s <- c(1L, 4L); e <- c(3L, 6L)
  expect_equal(seg_count(x, s, e), c(2, 2))
  expect_equal(seg_mean_na(x, s, e)$mean, c(2, 5))
  expect_equal(seg_mean_na(x, s, e)$n, c(2, 2))
})

test_that("seg_index and seg_id agree with the segment table", {
  s <- c(1L, 5L, 9L); e <- c(4L, 8L, 12L)
  expect_equal(seg_index(s, e), 1:12)
  expect_equal(seg_id(s, e), rep(1:3, each = 4))
})

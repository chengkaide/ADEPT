# Tests for the plateau filtering cascade.
#
# The direction filter is the one that used to be brittle: it demanded a
# perfectly monotonic sequence and collapsed to two plateaus when that failed.

seg_table <- function(ages, variance = NULL) {
  n <- length(ages)
  data.frame(
    Start = seq(0, by = 10, length.out = n),
    End   = seq(10, by = 10, length.out = n),
    Time_step = rep(10, n),
    Segment_Mean = ages,
    Variance = if (is.null(variance)) seq(0.01, 0.10, length.out = n) else variance,
    stringsAsFactors = FALSE
  )
}

test_that("longest_monotonic finds the longest increasing subsequence", {
  expect_equal(longest_monotonic(c(1, 2, 3), TRUE), 1:3)
  expect_equal(longest_monotonic(c(1, 3, 2, 4), TRUE), c(1L, 2L, 4L))
  expect_equal(longest_monotonic(c(3, 2, 1), TRUE), 1L)
  expect_equal(longest_monotonic(c(3, 2, 1), FALSE), 1:3)
  expect_equal(longest_monotonic(numeric(0)), integer(0))
  expect_equal(longest_monotonic(5), 1L)
})

test_that("longest_monotonic honours the tolerance", {
  a <- c(1, 2, 1.99, 3)             # the 1.99 is a 0.5 % reversal
  expect_equal(longest_monotonic(a, TRUE, tol = 0), c(1L, 2L, 4L))
  expect_equal(longest_monotonic(a, TRUE, tol = 0.02), 1:4)
})

test_that("a perfectly ascending sequence keeps every plateau", {
  s <- seg_table(c(21.6, 22.6, 23.6))
  out <- filter_direction(s, s$Segment_Mean, "Forward")
  expect_equal(out, c(21.6, 22.6, 23.6))
})

test_that("a sub-analytical reversal no longer collapses the result", {
  # A 0.5 % dip in the middle of an otherwise ascending stack. The old strict
  # test demanded perfect monotonicity, failed, and left only two plateaus.
  s <- seg_table(c(21.6, 22.0, 21.9, 23.6))
  out <- filter_direction(s, s$Segment_Mean, "Forward")
  expect_equal(sum(!is.na(out)), 4L)
  expect_equal(out[!is.na(out)], c(21.6, 22.0, 21.9, 23.6))
})

test_that("strict keeps reproducing the v1.1.0 behaviour", {
  s <- seg_table(c(21.6, 22.0, 21.9, 23.6))
  out <- filter_direction(s, s$Segment_Mean, "Forward", method = "strict")
  # not monotonic -> first plateau plus the smallest-variance remaining one
  expect_equal(sum(!is.na(out)), 2L)
  expect_equal(out[1], 21.6)
  # ... which is exactly the collapse v1.2.0 fixes
  expect_gt(sum(!is.na(filter_direction(s, s$Segment_Mean, "Forward"))), 2L)
})

test_that("a descending sequence is only kept in Reverse", {
  s <- seg_table(c(30, 25, 20))
  expect_equal(sum(!is.na(filter_direction(s, s$Segment_Mean, "Forward"))), 1L)
  expect_equal(sum(!is.na(filter_direction(s, s$Segment_Mean, "Reverse"))), 3L)
})

test_that("filter_age_range drops out-of-range and non-positive ages", {
  s <- seg_table(c(-5, 0, 100, 5000))
  expect_equal(is.na(filter_age_range(s, 0, 4540)), c(TRUE, FALSE, FALSE, TRUE))
})

test_that("filter_variance and filter_resolution gate on the right column", {
  s <- seg_table(c(1, 2, 3, 4), variance = c(0.01, 0.5, 0.01, 0.01))
  f1 <- filter_age_range(s, 0, 4540)
  expect_equal(is.na(filter_variance(s, f1, 0.1192)), c(FALSE, TRUE, FALSE, FALSE))

  s2 <- seg_table(c(1, 2, 3))
  s2$Time_step <- c(10, 2, 10)              # the first segment is always exempt
  f2 <- c(1, 2, 3)
  expect_equal(is.na(filter_resolution(s2, f2, min_res = 5)), c(FALSE, TRUE, FALSE))
})

test_that("apply_filters records the whole cascade and the counts", {
  s <- seg_table(c(21.6, 22.6, 23.6))
  out <- apply_filters(s, min_age = 0, max_age = 4540, var_threshold = 1,
                       min_res = 5, direction = "Forward")
  expect_true(all(c("Filter_1", "Filter_2", "Filter_3", "Filter_4",
                    "Final_Age", "Total_Segments", "Confirmed_Plateaus",
                    "Plateau_Numbers") %in% names(out)))
  expect_equal(out$Confirmed_Plateaus[1], 3L)
  expect_equal(out$Total_Segments[1], 3L)
  # the legacy name is kept as an alias
  expect_equal(out$Plateau_Numbers, out$Total_Segments)
  expect_equal(out$Final_Serial_Number, c(1L, 2L, 3L))
})

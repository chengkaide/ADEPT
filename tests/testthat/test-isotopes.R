# Regression tests for the built-in U-Pb age conversion.
#
# These pin down the decay constants and, more importantly, the *units* they
# are expressed in. An earlier revision divided a per-Ma uncertainty by a
# per-year decay constant and reported a relative uncertainty of 53500 %.

test_that("206Pb/238U ages match hand calculations", {
  # t = ln(1 + r) / lambda238
  expect_equal(age_206_238(0.1), log1p(0.1) / 1.55125e-10 / 1e6)
  expect_equal(age_206_238(0.5), 2613.796023, tolerance = 1e-5)
  expect_equal(age_206_238(0.00314), 20.210027, tolerance = 1e-5)
})

test_that("207Pb/235U ages match hand calculations", {
  expect_equal(age_207_235(0.6), 477.233720, tolerance = 1e-5)
  expect_equal(age_207_235(0.1), log1p(0.1) / 9.8485e-10 / 1e6)
})

test_that("207Pb/206Pb ages match published values", {
  expect_equal(age_207_206(0.05), 193.947601, tolerance = 1e-4)
  expect_equal(age_207_206(0.2), 2825.462288, tolerance = 1e-3)
  # below the present-day intercept there is no solution
  expect_true(is.na(age_207_206(0.01)))
  expect_true(is.na(age_206_238(-1)))
  expect_true(is.na(age_206_238(NA_real_)))
})

test_that("isotopes are consistent: a concordant point gives one age", {
  # 2 Ga: build the ratios from the forward equations, then invert
  t <- 2e9
  r68 <- exp(1.55125e-10 * t) - 1
  r75 <- exp(9.8485e-10 * t) - 1
  r76 <- (1 / 137.818) * (exp(9.8485e-10 * t) - 1) /
         (exp(1.55125e-10 * t) - 1)
  # the converters return Ma, so a 2 Ga point must come back as 2000
  expect_equal(age_206_238(r68), 2000, tolerance = 1e-6)
  expect_equal(age_207_235(r75), 2000, tolerance = 1e-6)
  expect_equal(age_207_206(r76), 2000, tolerance = 1e-4)
})

test_that("count and ratio inputs give the same ages", {
  a <- counts_to_ages(Pb206 = 100, Pb207 = 6, U238 = 1000)
  b <- ratios_to_ages(100 / 1000, 6 / (1000 / 137.818), 6 / 100)
  expect_equal(a$Age68, b$Age68)
  expect_equal(a$Age75, b$Age75)
  expect_equal(a$Age76, b$Age76)
})

test_that("decay-constant relative uncertainties are in the right units", {
  # A relative uncertainty is dimensionless: a value above a few percent is a
  # unit error, not a physical result.
  r68 <- adept_decay_rel_1s("Age68")
  r75 <- adept_decay_rel_1s("Age75")

  expect_equal(r68, 0.0005350, tolerance = 1e-3)
  expect_equal(r75, 0.0006803, tolerance = 1e-3)
  expect_true(r68 > 1e-5 && r68 < 1e-2)
  expect_true(r75 > 1e-5 && r75 < 1e-2)

  # 207Pb/206Pb is age dependent and less sensitive than the other two,
  # because the sensitivities to lambda235 and lambda238 partly cancel.
  r76_2ga <- adept_decay_rel_1s("Age76", 2000)
  expect_equal(r76_2ga, 0.001497, tolerance = 1e-3)
  expect_true(r76_2ga < r75 * 10)
  expect_false(isTRUE(all.equal(r76_2ga,
                                adept_decay_rel_1s("Age76", 100))))
})

test_that("adept_decay_rel_1s handles vectors and unknown systems", {
  v <- adept_decay_rel_1s(c("Age68", "Age75", "Age76"), c(100, 200, 2000))
  expect_length(v, 3)
  expect_true(all(is.finite(v)))
  expect_true(all(is.na(adept_decay_rel_1s("nonsense"))))
})

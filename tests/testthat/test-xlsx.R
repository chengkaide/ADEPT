# Tests for the built-in OOXML reader and writer.
#
# These are the parts that replace readxl / openxlsx / writexl, so a silent
# regression here would corrupt user data rather than merely annoy.

test_that("CRC-32 matches the standard test vectors", {
  expect_equal(crc32_raw(charToRaw("123456789")), 3421780262)
  expect_equal(crc32_raw(charToRaw("The quick brown fox jumps over the lazy dog")),
               1095738169)
  expect_equal(crc32_raw(raw(0)), 0)
})

test_that("column letters convert correctly", {
  expect_equal(col_letters_to_int("A"), 1L)
  expect_equal(col_letters_to_int("Z"), 26L)
  expect_equal(col_letters_to_int("AA"), 27L)
  expect_equal(col_letters_to_int("AZ"), 52L)
  expect_equal(col_letters_to_int("BA"), 53L)
  expect_equal(col_letters_to_int("ZZ"), 702L)
  expect_equal(col_letters_to_int("AAA"), 703L)
  expect_equal(.num_to_col_letters(c(1L, 26L, 27L, 702L, 703L)),
               c("A", "Z", "AA", "ZZ", "AAA"))
})

test_that("the round trip preserves numbers, NA, text and odd characters", {
  df1 <- data.frame(Analysis = c("Z-01", "Z-02"),
                    Age = c(1234.5678, NA),
                    Note = c("ok", "a <b> & 'c'"),
                    stringsAsFactors = FALSE)
  df1$Age <- as.numeric(df1$Age)
  df2 <- data.frame(x = 1:5,
                    y = c(0.1, 1e-10, 3.14159, NA, -2.5))

  f <- tempfile(fileext = ".xlsx")
  on.exit(unlink(f), add = TRUE)
  write_xlsx_base(list(Summary = df1, Full = df2), f, engine = "internal")

  back <- adept_read(f)
  expect_equal(names(back), c("Summary", "Full"))
  expect_equal(colnames(back$Summary), c("Analysis", "Age", "Note"))
  expect_equal(back$Summary$Analysis, c("Z-01", "Z-02"))
  expect_equal(back$Summary$Age, c(1234.5678, NA))
  expect_equal(back$Summary$Note, c("ok", "a <b> & 'c'"))
  expect_equal(back$Full$x, 1:5)
  expect_equal(back$Full$y, c(0.1, 1e-10, 3.14159, NA, -2.5))
})

test_that("written files are readable by readxl when it is installed", {
  skip_if_not_installed("readxl")
  df <- data.frame(a = c(1.5, 2.5), b = c("x", "y"), stringsAsFactors = FALSE)
  f <- tempfile(fileext = ".xlsx")
  on.exit(unlink(f), add = TRUE)
  write_xlsx_base(list(Sheet1 = df), f, engine = "internal")

  got <- as.data.frame(readxl::read_excel(f, sheet = "Sheet1"))
  expect_equal(colnames(got), c("a", "b"))
  expect_equal(got$a, c(1.5, 2.5))
  expect_equal(as.character(got$b), c("x", "y"))
})

test_that("the bundled example workbook reads back consistently", {
  f <- system.file("extdata", "Input(1sample).xlsx", package = "ADEPT")
  skip_if(f == "", "example workbook not installed")

  d <- adept_read(f)
  expect_equal(names(d), "Sheet1")
  expect_equal(nrow(d[[1]]), 410)
  expect_true(all(c("Analysis", "Time", "Pb206_U238", "Age68") %in%
                    colnames(d[[1]])))
  expect_equal(d[[1]]$Analysis[1], "DJ27_18")
  expect_true(is.numeric(d[[1]]$Time))
})

test_that("adept_read falls back to CSV and rejects legacy xls", {
  f <- tempfile(fileext = ".csv")
  on.exit(unlink(f), add = TRUE)
  write.csv(data.frame(a = 1:3, b = letters[1:3]), f, row.names = FALSE)
  d <- adept_read(f)
  expect_equal(nrow(d[[1]]), 3)

  expect_error(adept_read(tempfile(fileext = ".xls")), "not found")
})

test_that("xml helpers escape and unescape symmetrically", {
  s <- c("a & b", "x < y", "q \"q\"", "it's")
  expect_equal(xml_unescape(xml_escape(s)), s)
})

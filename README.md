# ADEPT: Automated Depth Profiling Technique

An R package for automated identification and extraction of age plateaus from
U-Pb depth profiling data.

**Reference:** Wang et al. (2026). Quantitative depth profiling of zircon
crystallization histories: New insights into Himalayan magmatic-tectonic
evolution. *JGR: Solid Earth, 131*, e2025JB033324.

## Installation

### From GitHub (recommended)

```r
if (!require("remotes")) install.packages("remotes")
remotes::install_github("yaowang-space/ADEPT")
```

### From tarball (no internet required after download)

1. Download `ADEPT_1.0.1.tar.gz` from [GitHub Releases](https://github.com/yaowang-space/ADEPT/releases)
2. Install locally:

```r
install.packages("ADEPT_1.0.1.tar.gz", repos = NULL, type = "source")
```

### From local source

```r
remotes::install_local("path/to/ADEPT")
```

### Dependencies

**ADEPT has no hard external dependencies.** Every part of the pipeline runs on
base R:

| Original dependency | Now handled by |
|---|---|
| `ggplot2` | base `graphics` / `grDevices` |
| `changepoint` | built-in PELT (`R/pelt.R`) |
| `forecast` | `stats::arima()` + KPSS/AIC order search |
| `IsoplotR` | built-in U-Pb decay equations (`R/isotopes.R`) |
| `readxl`, `openxlsx` | built-in OOXML reader (`R/xlsx.R`) |
| `writexl` | built-in OOXML writer, *or* `writexl` if installed |
| `zoo`, `dplyr`, `MASS`, `gridExtra` | never used — removed |

Optional packages, all in `Suggests`:

- `writexl` — 10-100x faster Excel writing (the built-in writer is used when it
  is absent; force it with `options(ADEPT.xlsx_engine = "internal")`)
- `mcp` + JAGS — Bayesian MCMC posterior analysis
- `shiny` (+ `DT`) — graphical interface

## Graphical interface

No R coding needed — launch the Shiny workbench:

```r
library(ADEPT)
adept_gui()
```

It lets you upload an `.xlsx` file, adjust every parameter, page through the
depth profiles, and download the Excel workbook or a combined PDF of all
profiles. The app source lives in
[`inst/shiny/app.R`](inst/shiny/app.R) — edit it directly, no reinstall
required. See [`inst/shiny/README_gui.md`](inst/shiny/README_gui.md) for a
walkthrough.

## Quick Start

```r
library(ADEPT)

# Load example data
example_file <- system.file("extdata", "Input(1sample).xlsx", package = "ADEPT")

# Default usage: Forward direction, no MCMC
result <- adept(example_file)

# Reverse direction (for descending age sequences), with MCMC
result <- adept(example_file, filter_direction = "Reverse", mcmc = TRUE)

# Access results
head(result$summary)   # Simplified summary (confirmed plateaus only)
head(result$full)      # Full results (all plateaus, all columns)
result$plots[[1]]      # Display first depth profile plot
```

### Example Output

![Example depth profile plot](man/figures/example_plot.png)

*Age vs. ablation time with LOESS-smoothed trend (black line) and identified
age plateaus (red segments with uncertainty bands).*

## Main Function

```r
adept(
  file_path,                   # Path to input Excel file
  chunk_size             = 411,
  lower_ablation_time    = 29,     # seconds
  upper_ablation_time    = 58,     # seconds
  max_age_limit          = 4540,   # Ma
  min_age_limit          = 0,      # Ma
  min_plateau_resolution = NULL,   # NULL = auto (>= 5 s)
  variance_threshold     = 0.1192,
  filter_direction       = "Forward",  # "Forward" or "Reverse"
  outlier_method         = "arima",    # "arima" | "mad" | "none"
  outlier_sd             = 2,          # residual SD threshold
  u238u235               = 137.818,    # for count-based input
  mcmc                   = FALSE,      # Bayesian MCMC?
  make_plots             = TRUE,       # Build depth profile plots?
  save_plots_to_disk     = TRUE,       # Write PDFs to plot_dir?
  output_path            = NULL,       # NULL = auto-named; NA = do not write
  plot_dir               = NULL,       # NULL = same as input
  keep_profiles          = FALSE,      # Return per-zircon data for re-plotting
  progress               = NULL,       # Callback fn(fraction, detail)
  verbose                = TRUE
)
```

> `plot` is still accepted as a deprecated alias for `make_plots`.

## Parameter Details

| Parameter | Type | Default | Description |
|---|---|---|---|
| `file_path` | character | (required) | Path to input `.xlsx` file |
| `chunk_size` | numeric | 411 | Rows per processing chunk |
| `lower_ablation_time` | numeric | 29 | Min ablation time in seconds |
| `upper_ablation_time` | numeric | 58 | Max ablation time in seconds |
| `max_age_limit` | numeric | 4540 | Maximum valid age (Ma) |
| `min_age_limit` | numeric | 0 | Minimum valid age (Ma) |
| `min_plateau_resolution` | numeric / NULL | NULL | Min plateau duration (s). NULL defaults to 5 |
| `variance_threshold` | numeric | 0.1192 | Maximum plateau variance |
| `filter_direction` | character | "Forward" | "Forward" = keep ascending; "Reverse" = keep descending |
| `direction_method` | character | "monotonic" | "monotonic" = longest monotonic subsequence (tolerates noise); "strict" = v1.1.0 behaviour |
| `direction_tolerance` | numeric | 0.02 | Relative reversal tolerated by `direction_method = "monotonic"` |
| `outlier_method` | character | "arima" | "arima" (published method), "mad" (fast), "none" |
| `outlier_sd` | numeric | 2 | Threshold in residual SDs |
| `preprocess` | character | "arima_loess" | "arima_loess" (published pipeline) or "robust_loess" |
| `smooth` | character | "loess" | `"none"` skips smoothing entirely - use it when the input is already corrected for down-hole fractionation |
| `calibration_uncertainty` | numeric | 0.03 | Relative 1-sigma reproducibility of the primary reference material |
| `u238u235` | numeric | 137.818 | 238U/235U for count-based input |
| `validate_input` | logical | `TRUE` | Run the input sanity checks and report actionable problems? |
| `mcmc` | logical | `FALSE` | Run MCMC posterior analysis? |
| `make_plots` | logical | `TRUE` | Build depth profile plots? |
| `save_plots_to_disk` | logical | `TRUE` | Write plot PDFs to `plot_dir`? |
| `output_path` | character / NULL / NA | `NULL` | Output Excel path. `NULL` = auto, `NA` = do not write |
| `plot_dir` | character / NULL | `NULL` | Plot PDF directory. `NULL` = directory of the input file |
| `keep_profiles` | logical | `FALSE` | Return per-zircon series + plateau tables as `$profiles` |
| `progress` | function / NULL | `NULL` | `function(fraction, detail)` progress callback |
| `verbose` | logical | `TRUE` | Print per-zircon messages? |

## Input Data Format

The input Excel file must contain sheets with at least one of the following:

### Format 1: Direct ages
Columns: `Analysis`, `Time`, `Age68`, `Age75`, `Age76`

### Format 2: Raw isotope counts
Columns: `Analysis`, `Time`, `Pb206`, `Pb207`, `U238`

### Format 3: Isotopic ratios
Columns: `Analysis`, `Time`, `Pb206_U238`, `Pb207_U235`, `Pb207_Pb206`

### Format 4: Direct ages with per-point uncertainties
Columns: `Analysis`, `Time`, `Age68`, plus `Age68_1s`
(and optionally `Age75_1s`, `Age76_1s`)

Adding per-point 1-sigma columns changes how a plateau is reduced:

- the **plateau age** becomes the inverse-variance weighted mean of the
  measured ages instead of the mean of the LOESS-smoothed curve. The weights
  already carry the information the smoothing was there to supply.
- the **plateau uncertainty** comes from the weights (`1/sqrt(sum(w))`) rather
  than from the scatter of the points.
- **MSWD** is reported, together with its probability. Around 1 means the
  spread agrees with the errors supplied. Much above 1 means either real age
  heterogeneity inside the plateau or sigmas that are too small - the usual
  cause is counting statistics quoted without the external reproducibility.
  Much below 1 means sigmas that are too large.
- `Points with uncertainty` reports how many points actually entered the
  calculation; points with a missing or non-positive sigma are dropped.

The suffix is **1-sigma, not 2-sigma**. Passing 2-sigma makes MSWD four times
smaller, which turns a failing plateau into an apparently acceptable one.

If those sigmas already include the external reproducibility of the reference
material, set `calibration_uncertainty = 0` so that term is not counted twice.
If the data were also corrected for down-hole fractionation before they
reached ADEPT, set `smooth = "none"` as well: such a profile is flat inside a
domain, and fitting LOESS to it a second time would round off the real domain
boundaries.

Any additional numeric columns (e.g., trace elements) will be automatically
detected and their plateau means will be included in the output.

### Pairing Format 4 with an external reduction

The combination this format was designed for:

```
raw counts
  --[external reduction: dead time, common Pb, F(tau), jackknife]-->
per-window ages + 1-sigma
  --[ADEPT, Format 4]-->
plateau ages + MSWD
```

Two settings matter when the data arrive that way:

- **`lower_ablation_time = 0`, with `upper_ablation_time` past the last
  window.** The reduction has already picked its own ablation interval, and
  ADEPT's default 29-58 s window would cut a piece off the end. The window is
  deliberately *not* skipped automatically: a per-point sigma column says
  nothing about where the window ends, and guessing wrong silently changes how
  many points enter the segmentation.
- **`calibration_uncertainty = 0`** if those sigmas already include the
  external reproducibility of the reference material. Most standard-based
  reductions do include it. Counting it twice roughly multiplies the reported
  uncertainty by ten on a 460 Ma zircon (13.7 Ma instead of 1.8 Ma), which
  buries everything the per-point errors were there to resolve.

`smooth = "none"` is worth trying too, since it segments the measured profile
rather than a smoothed version of it. On data that are already window-averaged
the difference is usually small - the windows are smooth to begin with - but it
removes a step that is not needed.

A sheet may also stack as many zircons as you like: rows are grouped by
`Analysis`, and each group is segmented on its own. `Time` only has to increase
*within* a group, so every zircon's clock can restart from zero.

## Output

### Excel file (two sheets)

**Sheet 1 — Summary:** Confirmed age plateaus only.

| Column | Description |
|---|---|
| Analysis | Sample analysis ID |
| Group | Chunk group number |
| Points | Number of data points in chunk |
| Final serial number | Plateau rank (1, 2, ...) |
| Integration time | Plateau duration (s) |
| Final age (Ma) | Plateau age |
| Final total uncertainty (Ma) | Combined uncertainty |
| Concordance (%) | Age68/Age75 × 100 |
| Pb206/U238 age mean (Ma) | Mean 206Pb/238U age |
| Pb206/U238 total uncertainty (Ma) | Its uncertainty |
| Pb207/U235 age mean (Ma) | Mean 207Pb/235U age |
| Pb207/U235 total uncertainty (Ma) | Its uncertainty |
| Pb207/Pb206 age mean (Ma) | Mean 207Pb/206Pb age |
| Pb207/Pb206 total uncertainty (Ma) | Its uncertainty |
| `*_Mean` | Means of extra input columns |

**Sheet 2 — Full Results:** All plateaus with complete statistics including
MCMC posterior estimates (if `mcmc = TRUE`), slope/intercept, and filter flags.

```r
list(
  summary  = data.frame,   # Sheet 1 content
  full     = data.frame,   # Sheet 2 content
  plots    = list(),       # ggplot objects (if make_plots = TRUE)
  profiles = list()        # per-zircon data + segments (if keep_profiles = TRUE)
)
```

All numeric columns are returned as `numeric` (and written to Excel as
numbers, not text).

### Return Value

## Processing Pipeline

```
Input Excel → Format Detection → ARIMA Outliers → Discordance Filter
→ Mean Fill → LOESS Smoothing → Standardization → PELT Segmentation
→ Plateau Statistics → 4-Step Filtering → [MCMC] → Output
```

## Development

The package ships a regression test suite and a CI workflow. Both were added
because every bug found in this code so far was found by hand: a PELT pruning
error at `minseglen > 1`, and a per-Ma/per-year mix-up in the decay-constant
uncertainty.

```r
# from a source checkout
testthat::test_local()
```

`tests/testthat/` covers the parts that are easy to get subtly wrong and hard
to notice:

| file | what it pins down |
|---|---|
| `test-isotopes.R` | the U-Pb conversions against hand calculations, and that the decay-constant uncertainties are relative (dimensionless) values |
| `test-pelt.R` | changepoint vectors confirmed against `changepoint::cpt.mean()`, and that every segment respects `minseglen` |
| `test-segstats.R` | the vectorised prefix-sum statistics against the `sapply` + `lm()` loops they replaced |
| `test-filtering.R` | the four-step cascade, including the two-plateau collapse that v1.2.0 removed |
| `test-xlsx.R` | the base-R OOXML reader and writer round trip |
| `test-adept.R` | end-to-end runs on the bundled example workbook, with the published ages hard-coded as the anchor |

`.github/workflows/R-CMD-check.yaml` runs `R CMD check` on Ubuntu, macOS and
Windows, for the current and the previous R release.

## Changelog

### Unreleased

**`smooth = "none"`** (new). Skips the smoothing step entirely and segments the
profile as measured. Use it when the input has already been corrected for
down-hole fractionation - an F(tau) correction applied by an external
reduction, say. Such a profile is flat inside a domain, so fitting LOESS to it
again would round off the real domain boundaries. The outlier screen selected
by `preprocess` still runs, and the profile is still scaled before PELT; only
the smoothing is skipped. The default stays `"loess"`, so existing results are
unchanged - the bundled example still gives 21.61724 / 22.56383 / 23.60464 Ma.

Also fills in five parameters that v1.2.0 introduced but the parameter table
below never listed: `preprocess`, `calibration_uncertainty`,
`direction_method`, `direction_tolerance` and `validate_input`.

**Per-point uncertainties, and MSWD** (new, Format 4). Supplying `Age68_1s`
(or `Age75_1s` / `Age76_1s`) alongside the ages switches the plateau age to an
inverse-variance weighted mean of the measured ages, and adds `MSWD`,
`MSWD probability` and `Points with uncertainty` to the output. MSWD is the
number the U-Pb community uses to decide whether a plateau is internally
consistent: around 1 means the scatter matches the quoted errors.

This also closes a gap the package had from the start. Until now `Plateau
uncertainty` was the scatter of the points divided by sqrt(n), which measures
precision but says nothing about whether the points agree with one another -
a plateau built from two discordant populations could look precise. MSWD is
the test for that, and it needs per-point errors to exist.

Without sigma columns nothing changes: the three new columns are present and
`NA`, and every published number is reproduced exactly.

### 1.2.0

Focused on the two things a reviewer is most likely to challenge: where the
reported uncertainty comes from, and how much of the result rests on the
parameter choices. Existing numbers do not move - every change is either
additive or opt-in.

**Uncertainty budget.** The reference-material reproducibility was a literal
`0.03` buried in `calc_uncertainty()`, and on the bundled example it accounted
for **93-95 % of the reported total uncertainty** (the measured within-plateau
scatter contributed only 5-7 %). It is now the parameter
`calibration_uncertainty`, and the budget is decomposed:

| column | meaning |
|---|---|
| `Random uncertainty (Ma)` | within-plateau scatter + reference-material reproducibility |
| `Systematic uncertainty (Ma)` | decay constants (fully correlated between samples) |
| `Decay constant uncertainty (Ma)` | the decay-constant term on its own |
| `Total uncertainty incl. decay (Ma)` | quadrature sum of the two |
| `Decay system` | which system the plateau age came from |

The legacy `Total uncertainty` deliberately keeps excluding the decay
constants, so previously published numbers are unchanged. Decay-constant
relative 1-sigma values follow IsoplotR (Jaffey et al. 1971; Hiess et al.
2012): 0.0535 % for 206Pb/238U, 0.0680 % for 207Pb/235U, and ~0.15 % at 2 Ga
for 207Pb/206Pb (age-dependent, evaluated per plateau - the two decay constants
partly cancel, which is why 207Pb/206Pb ages are less sensitive).

**Parameter sensitivity.** New `adept_sensitivity()` sweeps
`variance_threshold`, `min_plateau_resolution` and `filter_direction`, and
returns a tidy table of plateau counts and ages. On the bundled example the
plateau count runs 0 / 1 / 3 for thresholds 0.05 / 0.08 / 0.1192 and 3 / 1 / 0
for minimum durations of 5 / 8 / 15 s - i.e. the defaults sit on a cliff, which
is exactly what a supplementary table should document.

**Plateau selection.** `filter_direction()` used `all(diff(ages) >= 0)`, a
perfect-monotonicity test that fails on 11.5 % of 3-plateau and 26.5 % of
5-plateau sequences at only 2 % noise, and then collapsed to "first plateau +
smallest variance" (at most two). It now keeps the **longest monotonic
subsequence** with a tolerance (`direction_method = "monotonic"`, default), so
every plateau consistent with the trend is retained.
`direction_method = "strict"` restores the old behaviour exactly.

**Input validation.** `validate_segment_data()` reports missing columns, a Time
axis that looks like milliseconds rather than seconds, non-monotonic time, and
negative or older-than-Earth ages *inside the ablation window*. Structural
problems skip the zircon with an actionable message; value problems warn only.

**Robust pre-processing (opt-in).** `preprocess = "robust_loess"` drops the
ARIMA outlier screen and fits `loess(family = "symmetric")`, so outliers are
downweighted by an M-estimator instead of being deleted. Removes the least
stable step in the pipeline.

**Other.** `Plateau_Numbers` is now also available as `Total_Segments` (the old
name was misleading: it is the number of PELT segments, not of confirmed
plateaus); `Confirmed plateaus` added. Workaround note: `parallel` added to
`Suggests` for multi-core sweeps.

### 1.1.0

**Dependency removal.** The package now runs on base R alone. `ggplot2`,
`changepoint`, `forecast`, `IsoplotR`, `readxl` and `openxlsx` are gone;
`writexl` became optional. This insulates the results from upstream changes —
verified against the previous implementation, segment boundaries and ages are
identical (differences ≤ 4e-12, i.e. floating-point noise).

- Built-in **OOXML reader/writer** (`R/xlsx.R`). Reading is validated against
  `readxl` on the bundled example (identical column names and values, max
  difference 0) and is ~1.7x slower than `readxl`; writing produces files that
  Excel, LibreOffice and `readxl` all read back correctly.
- Built-in **PELT** changepoint detection (`R/pelt.R`), L2 cost. Reproduces
  `changepoint::cpt.mean(penalty = "Manual", minseglen = 1)` changepoint sets
  **exactly** on all tested series. The AIC penalty (always 4 for a univariate
  mean model) is now a package constant instead of a package lookup.
- Built-in **U-Pb age conversion** (`R/isotopes.R`), matching `IsoplotR` to
  < 1e-4 Ma over the whole geological range.
- Built-in **ARIMA order selection** (`stats::arima()` + KPSS + AIC grid),
  replacing `forecast::auto.arima()`. New argument `outlier_method =
  c("arima", "mad", "none")`.
- **Vectorised segment statistics** (`R/segstats.R`): prefix sums instead of
  `sapply` + `lm()` per segment. Segment means 71x faster, segment regression
  ~2000x faster, and the whole pipeline ~2x faster end-to-end on a 40-zircon
  workbook — while dropping a C dependency.
- **Base-graphics plotting**: `plot_depth_profile()` now returns an
  `adept_profile` record drawn by `draw_profile()`.

### 1.0.2

- **Fix:** numeric results were written to Excel as *text*. `adept()` built
  each row with `c()`, which coerces numbers to character. Rows are now named
  lists, so every numeric column keeps its type.
- **Fix:** `mcmc = TRUE` without the `mcp` package crashed with
  `length of 'dimnames' [2] not equal to array extent`. MCMC columns are now
  always present (filled with `NA` when unavailable).
- **Fix:** `mcp_step` logic (`sum(!is.na(mcp_step)) >= 3` was a tautology) and
  the single-plateau case, which errored on a zero-length replacement.
- **Fix:** renamed `Final Serial Number` → `Final serial number` so the column
  name no longer changes with `mcmc`.
- **Fix:** `tools::file_path_sans_ext` used without declaring `tools` in
  `Imports`; `IsoplotR` is now checked with an actionable error message.
- **Robustness:** `auto.arima`, LOESS and PELT are wrapped in guards, so a
  single bad zircon no longer aborts the whole run — it emits an empty row.
- **Removed** unused dependencies `dplyr`, `MASS`, `gridExtra`, `zoo`.
- **Added** `man/` documentation, a Shiny GUI (`adept_gui()`), the
  `keep_profiles` and `progress` arguments, and `make_plots` /
  `save_plots_to_disk` (replacing the ambiguous `plot`).

## Citation

Wang, Y., Liu, Z., Zhang, L.-L., Liu, L., Walters, J., Cawood, P. A.,
Yang, L., Wang, Q., Stockli, D. F., & Zhu, D.-C. (2026). Quantitative
depth profiling of zircon crystallization histories: New insights into
Himalayan magmatic-tectonic evolution. *Journal of Geophysical Research:
Solid Earth*, *131*, e2025JB033324. https://doi.org/10.1029/2025JB033324

## License

MIT

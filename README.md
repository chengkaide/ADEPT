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
| `outlier_method` | character | "arima" | "arima" (published method), "mad" (fast), "none" |
| `outlier_sd` | numeric | 2 | Threshold in residual SDs |
| `u238u235` | numeric | 137.818 | 238U/235U for count-based input |
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
Columns: `Analysis`, `Time`, `Pb206`, `Pb207`, `U235`, `U238`

### Format 3: Isotopic ratios
Columns: `Analysis`, `Time`, `Pb206_U238`, `Pb207_U235`, `Pb207_Pb206`

Any additional numeric columns (e.g., trace elements) will be automatically
detected and their plateau means will be included in the output.

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

## Changelog

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

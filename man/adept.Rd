\name{adept}
\alias{adept}
\title{Automated Depth Profiling Technique}
\usage{
adept(
  file_path,
  chunk_size = 411,
  lower_ablation_time = 29,
  upper_ablation_time = 58,
  max_age_limit = 4540,
  min_age_limit = 0,
  min_plateau_resolution = NULL,
  variance_threshold = 0.1192,
  filter_direction = c("Forward", "Reverse"),
  direction_method = c("monotonic", "strict"),
  direction_tolerance = 0.02,
  outlier_method = c("arima", "mad", "none"),
  outlier_sd = 2,
  preprocess = c("arima_loess", "robust_loess"),
  calibration_uncertainty = 0.03,
  u238u235 = ADEPT_U238U235,
  validate_input = TRUE,
  mcmc = FALSE,
  make_plots = TRUE,
  save_plots_to_disk = TRUE,
  output_path = NULL,
  plot_dir = NULL,
  keep_profiles = FALSE,
  progress = NULL,
  verbose = TRUE,
  plot = TRUE
)
}
\arguments{
\item{file_path}{Path to the input Excel (.xlsx) or CSV file.}

\item{chunk_size}{Number of rows per processing chunk (default 411).}

\item{lower_ablation_time}{Minimum effective ablation time in seconds (default 29).}

\item{upper_ablation_time}{Maximum effective ablation time in seconds (default 58).}

\item{max_age_limit}{Maximum valid age in Ma (default 4540).}

\item{min_age_limit}{Minimum valid age in Ma (default 0).}

\item{min_plateau_resolution}{Minimum plateau duration in seconds.
\code{NULL} or values < 5 default to 5 seconds. The first segment is exempt,
because the onset of the ablation signal is usually truncated by
\code{lower_ablation_time}.}

\item{variance_threshold}{Maximum allowed intra-plateau variance (default 0.1192).}

\item{filter_direction}{\code{"Forward"} keeps ascending age sequences;
\code{"Reverse"} keeps descending age sequences.}

\item{direction_method}{\code{"monotonic"} (default) keeps the longest set of
plateaus whose ages follow the requested trend, tolerating reversals smaller
than \code{direction_tolerance}. \code{"strict"} reproduces the v1.1.0
behaviour, which required a perfectly monotonic sequence and otherwise fell
back to the first plateau plus the one with the smallest variance (at most two
plateaus).}

\item{direction_tolerance}{Relative reversal tolerated by
\code{direction_method = "monotonic"} (default 0.02, i.e. 2 \%).}

\item{outlier_method}{Outlier detection applied before smoothing:
\code{"arima"} (default, matches the published method), \code{"mad"}
(fast robust z-score) or \code{"none"}. Ignored when
\code{preprocess = "robust_loess"}.}

\item{outlier_sd}{Threshold, in residual standard deviations, beyond which a
point is discarded (default 2).}

\item{preprocess}{\code{"arima_loess"} (default, the published pipeline:
ARIMA residual screening followed by an ordinary LOESS fit) or
\code{"robust_loess"}, which skips the ARIMA screen entirely and fits
\code{loess(family = "symmetric")} so that outliers are downweighted rather
than deleted.}

\item{calibration_uncertainty}{Relative 1-sigma reproducibility of the primary
reference material (default 0.03, i.e. 3 \%). It dominated the reported total
uncertainty, so it is now an explicit parameter rather than a literal buried
in the code.}

\item{validate_input}{Logical. Run the input sanity checks (required columns,
monotonic and plausibly scaled Time axis, usable age values) and report
actionable problems? (default \code{TRUE}).}

\item{u238u235}{238U/235U ratio used to derive 235U for count-based input
(default 137.818).}

\item{mcmc}{Logical. Run Bayesian MCMC posterior analysis? (default \code{FALSE}).
Requires the suggested package \pkg{mcp} (and JAGS); if unavailable the MCMC
columns are \code{NA} rather than an error.}

\item{make_plots}{Logical. Build depth profile plots? (default \code{TRUE}).}

\item{save_plots_to_disk}{Logical. Write plot PDFs to \code{plot_dir}
(default \code{TRUE}); requires \code{make_plots = TRUE}.}

\item{output_path}{Output Excel path. \code{NULL} auto-generates from the input
name; \code{NA} skips writing the workbook.}

\item{plot_dir}{Directory for plot PDFs. \code{NULL} uses the directory of
\code{file_path}.}

\item{keep_profiles}{Logical. Also return the per-zircon smoothed series and
plateau tables under \code{$profiles}? (default \code{FALSE}).}

\item{progress}{Optional callback \code{function(fraction, detail)} invoked as
the run advances; used by \code{\link{adept_gui}}.}

\item{verbose}{Logical. Emit per-zircon progress messages? (default \code{TRUE}).}

\item{plot}{Deprecated alias for \code{make_plots}.}
}
\value{
Invisibly returns a list with \code{summary}, \code{full}, \code{plots} and
(optionally) \code{profiles}. Unless disabled, also writes a two-sheet Excel
file and PDF plots to disk.

All numeric columns are returned as numeric (earlier versions returned
everything as character, which made the exported workbook unusable for
further calculation). Plot objects are \code{adept_profile} records that are
drawn by the package's base-graphics renderer, so \pkg{ggplot2} is not needed.
}
\description{
Quantitatively identifies and extracts age plateaus from LA-ICP-MS U-Pb depth
profiling data. Supports three input formats: direct ages
(Age68/Age75/Age76), raw isotope counts (Pb206/Pb207/U238), or isotopic ratios
(Pb206_U238/Pb207_U235/Pb207_Pb206).
}
\details{
All core functionality runs on base R. OOXML reading and writing, U-Pb age
conversion, PELT changepoint segmentation, ARIMA order selection and plotting
are implemented inside the package, so results do not drift with changes in
third-party packages. \pkg{writexl}, \pkg{mcp} and \pkg{shiny} are optional.
}
\examples{
\dontrun{
result <- adept("Input.xlsx")
result <- adept("Input.xlsx", filter_direction = "Reverse")
result <- adept("Input.xlsx", mcmc = TRUE)
}
}

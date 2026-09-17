\name{adept_sensitivity}
\alias{adept_sensitivity}
\title{Sweep ADEPT parameters and tabulate the outcome}
\usage{
adept_sensitivity(
  file_path,
  variance_threshold = c(0.05, 0.08, 0.1192, 0.15, 0.25),
  min_plateau_resolution = c(5, 8, 10, 15),
  filter_direction = c("Forward", "Reverse"),
  calibration_uncertainty = 0.03,
  direction_method = c("monotonic", "strict"),
  direction_tolerance = 0.02,
  n_cores = 1L,
  progress = NULL,
  verbose = TRUE,
  ...
)
}
\arguments{
\item{file_path}{Path to the input workbook.}

\item{variance_threshold}{Numeric vector of variance thresholds to try.}

\item{min_plateau_resolution}{Numeric vector of minimum plateau durations
(seconds) to try.}

\item{filter_direction}{Character vector of directions to try.}

\item{calibration_uncertainty}{Relative 1-sigma reference-material
reproducibility, passed through to \code{\link{adept}}.}

\item{direction_method, direction_tolerance}{Passed through to
\code{\link{adept}}.}

\item{n_cores}{Number of workers. 1 (default) keeps everything in the current
process; > 1 uses \pkg{parallel}, which ships with R.}

\item{progress}{Optional callback \code{function(fraction, detail)}.}

\item{verbose}{Logical. Print progress?}

\item{...}{Further arguments passed to \code{\link{adept}}.}
}
\value{
A data.frame with one row per combination and columns
\code{variance_threshold}, \code{min_plateau_resolution},
\code{filter_direction}, \code{n_plateaus}, \code{ages} (semicolon separated)
and \code{mean_age}.
}
\description{
Runs \code{\link{adept}} once per parameter combination and returns a tidy
table of the plateau count and ages. Nothing is written to disk.
}
\details{
The plateau count and the resulting ages are sensitive to
\code{variance_threshold} and \code{min_plateau_resolution}. On the bundled
example the plateau count runs 0, 1, 3 for thresholds of 0.05, 0.08 and 0.1192,
and 3, 1, 0 for minimum durations of 5, 8 and 15 seconds. Reviewers routinely
ask how much of the interpretation rests on those choices, so this function
produces the table directly.
}
\examples{
\dontrun{
s <- adept_sensitivity("Input.xlsx")
print(s)

# parallel sweep on four cores
s <- adept_sensitivity("Input.xlsx", n_cores = 4)
}
}
\seealso{
\code{\link{adept}}
}

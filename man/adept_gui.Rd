\name{adept_gui}
\alias{adept_gui}
\title{Launch the ADEPT Shiny workbench}
\usage{
adept_gui(host = "127.0.0.1", port = NULL, launch.browser = TRUE, ...)
}
\arguments{
\item{host}{Host passed to \code{shiny::runApp} (default \code{"127.0.0.1"}).}

\item{port}{Port passed to \code{shiny::runApp}. \code{NULL} picks a free port.}

\item{launch.browser}{Open the system browser? (default \code{TRUE}).}

\item{...}{Further arguments forwarded to \code{shiny::runApp}.}
}
\value{
Invisibly returns the value of \code{shiny::runApp} (called for
its side effect). The Shiny app requires the suggested packages
\pkg{shiny} (mandatory) and \pkg{DT} (optional, for interactive tables).
}
\description{
Opens an interactive front-end for \code{\link{adept}}: upload an Excel
file, tweak every parameter, inspect the resulting depth profiles, and
download the Excel / PDF output.
}
\examples{
\dontrun{
if (requireNamespace("shiny", quietly = TRUE)) {
  adept_gui()
}
}
}

# ADEPT: Automated Depth Profiling Technique
# gui.R — Shiny front-end launcher
#
# The actual application lives in inst/shiny/app.R so that users can edit it
# directly without rebuilding the package. This file only exposes a launcher.

#' Launch the ADEPT Shiny workbench
#'
#' Opens an interactive front-end for \code{\link{adept}}: upload an Excel
#' file, tweak every parameter, inspect the resulting depth profiles, and
#' download the Excel / PDF output.
#'
#' @param host Host passed to \code{shiny::runApp} (default \code{"127.0.0.1"}).
#' @param port Port passed to \code{shiny::runApp}. \code{NULL} picks a free port.
#' @param launch.browser Open the system browser? (default \code{TRUE}).
#' @param ... Further arguments forwarded to \code{shiny::runApp}.
#'
#' @return Invisibly returns the value of \code{shiny::runApp} (called for
#'   its side effect). The Shiny app requires the suggested packages
#'   \pkg{shiny} (mandatory) and \pkg{DT} (optional, for interactive tables).
#'
#' @seealso \code{\link{adept}}
#' @export
#'
#' @examples
#' \dontrun{
#' if (requireNamespace("shiny", quietly = TRUE)) {
#'   adept_gui()
#' }
#' }
adept_gui <- function(host = "127.0.0.1", port = NULL,
                      launch.browser = TRUE, ...) {
  if (!requireNamespace("shiny", quietly = TRUE)) {
    stop("The ADEPT GUI needs the 'shiny' package.\n",
         "  install.packages(c(\"shiny\", \"DT\"))", call. = FALSE)
  }

  app_dir <- system.file("shiny", package = "ADEPT")
  if (!dir.exists(app_dir)) {
    stop("Could not locate the bundled Shiny app at ", app_dir,
         ". Reinstall ADEPT.", call. = FALSE)
  }

  if (!requireNamespace("DT", quietly = TRUE)) {
    message("Note: installing DT adds sortable, exportable result tables.\n",
            "  install.packages(\"DT\")")
  }

  shiny::runApp(appDir = app_dir, host = host, port = port,
                launch.browser = launch.browser, ...)
}

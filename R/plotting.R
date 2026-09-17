# ADEPT: Automated Depth Profiling Technique
# plotting.R — publication-ready depth profile plots in base graphics
#
# Replaces the ggplot2 dependency (and its ~10 transitive packages). The
# layout mirrors the original figures: raw ages as translucent dots, the
# LOESS trend as a black line, and the retained plateaus as red bars with
# uncertainty bands.

#' Pack a depth profile into a drawable object
#'
#' Returns a lightweight description of the plot. Drawing is deferred to
#' `draw_profile()` so that the same object can be sent to a PDF device, a
#' PNG device, or a Shiny `renderPlot()`.
#'
#' @param df data.frame with Time, Raw_Age, loess_Age columns
#' @param segments Segment data.frame with Start, End, Segment_Mean,
#'   Total_uncertainty and Filter_4
#' @param title Plot title (defaults to the Analysis name)
#' @return An object of class \code{adept_profile}
#' @keywords internal
plot_depth_profile <- function(df, segments, title = NULL,
                               label = NULL, title_text = NULL) {
  if (is.null(title) && !is.null(title_text)) title <- title_text
  if (is.null(title)) title <- paste("Analysis:", df[1, "Analysis"])
  structure(
    list(data = df, segments = segments, title = title, label = label),
    class = "adept_profile"
  )
}

#' Draw a depth profile on the current graphics device
#'
#' @param x An \code{adept_profile} object.
#' @param xlab,ylab Axis labels.
#' @param cex_point Point size for the raw ages.
#' @param ... Further arguments passed to \code{plot()}.
#' @return Invisibly returns `x`.
#' @keywords internal
draw_profile <- function(x, xlab = "Time (s)", ylab = "Age (Ma)",
                         cex_point = 0.55, ...) {
  df  <- x$data
  seg <- x$segments

  op <- graphics::par(mar = c(3.6, 3.8, 2.2, 1.0),
                      mgp = c(2.2, 0.55, 0),
                      tcl = -0.28,
                      xaxs = "r", yaxs = "r",
                      cex.axis = 0.85, cex.lab = 0.9, cex.main = 0.95)
  on.exit(graphics::par(op), add = TRUE)

  graphics::plot(df$Time, df$Raw_Age, type = "n",
                 xlab = xlab, ylab = ylab, main = x$title, ...)

  keep <- !is.na(seg$Filter_4)
  if (any(keep)) {
    lo <- seg$Segment_Mean[keep] - seg$Total_uncertainty[keep]
    hi <- seg$Segment_Mean[keep] + seg$Total_uncertainty[keep]
    graphics::rect(seg$Start[keep], lo, seg$End[keep], hi,
                   col = grDevices::adjustcolor("red", 0.15), border = NA)
  }

  graphics::points(df$Time, df$Raw_Age, pch = 16, cex = cex_point,
                   col = grDevices::adjustcolor("#f46f20", 0.35))
  graphics::lines(df$Time, df$loess_Age, lwd = 1.6)

  if (any(keep)) {
    graphics::segments(seg$Start[keep], seg$Segment_Mean[keep],
                       seg$End[keep],   seg$Segment_Mean[keep],
                       col = "red", lwd = 2)
  }
  graphics::box()
  invisible(x)
}

#' @export
print.adept_profile <- function(x, ...) {
  draw_profile(x, ...)
  invisible(x)
}

#' Save depth profile plots to PDF
#'
#' @param plots List of \code{adept_profile} objects.
#' @param output_dir Directory to save PDFs.
#' @param prefix Filename prefix (e.g. "Plot").
#' @param width_cm,height_cm Page size in cm (defaults match the original
#'   66 x 43.5 mm figure size).
#' @keywords internal
save_plots <- function(plots, output_dir = ".", prefix = "Plot",
                       width_cm = 6.6, height_cm = 4.35) {
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  for (i in seq_along(plots)) {
    lbl <- plots[[i]]$label
    if (is.null(lbl) || is.na(lbl) || !nzchar(lbl)) lbl <- i
    fname <- file.path(output_dir,
                       paste0(prefix, "_", make.names(lbl), ".pdf"))
    grDevices::pdf(fname, width = width_cm / 2.54, height = height_cm / 2.54)
    draw_profile(plots[[i]])
    grDevices::dev.off()
  }
  invisible(plots)
}

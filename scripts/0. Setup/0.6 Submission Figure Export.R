# 0.6 Submission Figure Export
# ---------------------------------------------------------------------------
# Writes journal submission ready copies of the manuscript figures, high
# resolution, etc.
#
# Requires the ragg package (grDevices::tiff needs a working cairo build,
# which is not available on every machine).
# ---------------------------------------------------------------------------

EXPORT_SUBMISSION_FIGURES <- TRUE            # set FALSE to skip the TIFF export
SUBMISSION_FIG_DIR        <- "figures/submission"
SUBMISSION_FIG_AUTHOR     <- "Fawcett"
SUBMISSION_PRINT_WIDTH_IN <- 6.5             # assumed final printed width (full text width)
SUBMISSION_TARGET_DPI     <- 800             # resolution at the printed width (journal: 800 preferred, 600 minimum)

save_submission_figure <- function(plot, number, width, height, units = "in",
                                   scale = 1,
                                   print_width = SUBMISSION_PRINT_WIDTH_IN,
                                   target_dpi  = SUBMISSION_TARGET_DPI,
                                   dir         = SUBMISSION_FIG_DIR,
                                   author      = SUBMISSION_FIG_AUTHOR,
                                   compression = "lzw") {
  if (!isTRUE(EXPORT_SUBMISSION_FIGURES)) return(invisible(NULL))
  if (!requireNamespace("ragg", quietly = TRUE)) {
    warning("Package 'ragg' is not installed; submission TIFF for Figure ", number, " not written.")
    return(invisible(NULL))
  }
  if (units != "in") stop("save_submission_figure() expects width and height in inches.")
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)

  drawn_w <- width * scale                     # physical size of the rendered canvas
  drawn_h <- height * scale
  dpi     <- max(300, ceiling(target_dpi * print_width / drawn_w))
  file    <- file.path(dir, sprintf("%s_Figure_%s.tif", author, number))

  ggplot2::ggsave(file, plot, device = ragg::agg_tiff,
                  width = width, height = height, units = "in", scale = scale,
                  dpi = dpi, compression = compression, bg = "white")

  px_w <- round(drawn_w * dpi); px_h <- round(drawn_h * dpi)
  message(sprintf("Submission figure written: %s (%d x %d px; %d dpi at %.1f in)",
                  file, px_w, px_h, round(px_w / print_width), print_width))
  invisible(file)
}

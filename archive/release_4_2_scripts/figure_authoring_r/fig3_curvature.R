#!/usr/bin/env Rscript

script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1]])
source(file.path(dirname(normalizePath(script_file)), "_helpers.R"))

cfg <- list(
  name = "fig3_curvature",
  pdf_width_pt = 504,
  pdf_height_pt = 346.08,
  png_window_width = 1400,
  png_window_height = 960,
  png_scale = 1,
  png_use_wrapper = TRUE
)

edits <- list(
  text_edit("Figure 3. Discrete curvature", x = 28, y = 29, size = 17),
  text_edit("All listed allowed values", x = 28, y = 48, size = 12.5),
  text_edit("curvature ΩK", x = 447, y = 293, size = 13),
  # Keep deletion labels visibly inside the plotting rectangle, off its rule.
  text_edit("N = 4 deleted", x = 439, y = 252, size = 12.5, anchor = "end"),
  text_edit("N = 6 deleted", x = 536, y = 252, size = 12.5, anchor = "end"),
  text_edit("2025 set", x = 28, y = 314, size = 12.5),
  text_edit("Nearest allowed value", x = 28, y = 336, size = 12.5),
  text_edit("The 2025 χ² minimum", x = 28, y = 358, size = 12.5),
  text_edit("No compact reading", x = 28, y = 393, size = 12.5),
  # Legend collided with the row labels: "DESI: +0.0023 +/- 0.0011" spans
  # x 66..206 at y 159 and "2024 S3" spans x 196..236 at y 150. Ten px of
  # horizontal overlap, nine px apart vertically. There is a 69px empty band
  # between the subtitle at y 48 and the plot frame top at y 94, and the
  # legend sits left of the frame (x < 250), so lift the whole block 51px.
  text_edit("data bands", y = 66, exact = TRUE),
  text_edit("Planck: \u22120.0106 \u00b1 0.0065", y = 88, exact = TRUE),
  text_edit("DESI: +0.0023 \u00b1 0.0011", y = 108, exact = TRUE)
)

run_figure(cfg, edits)

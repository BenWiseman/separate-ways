#!/usr/bin/env Rscript

script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1]])
source(file.path(dirname(normalizePath(script_file)), "_helpers.R"))

cfg <- list(
  name = "fig4_tilt",
  pdf_width_pt = 504,
  pdf_height_pt = 382.08,
  png_window_width = 1400,
  png_window_height = 1060,
  png_scale = 1,
  png_use_wrapper = TRUE
)

edits <- list(
  text_edit("Figure 4. Primordial tilt", x = 28, y = 29, size = 17),
  text_edit("All estimates are base", x = 28, y = 48, size = 12.5),
  text_edit("solid: dimension-zero", x = 66, y = 74, size = 12.5),
  text_edit("dashed: amplitude branch", x = 66, y = 90, size = 12.5),
  # Recover a small left gutter without reducing the longest data label.
  text_edit("Planck 2018 TT,TE,EE+lowE+lensing", x = 212, size = 12.5),
  # Separate the axis title from the explanatory notes at final page scale.
  text_edit("scalar spectral index", x = 355, y = 416, size = 13),
  text_edit("Prediction: nₛ", x = 28, y = 442, size = 12.5),
  text_edit("The two values are not", x = 28, y = 465, size = 12.5),
  text_edit("not a formal exclusion", x = 28, y = 483, size = 12.5),
  text_edit("The CPT-symmetric structure", x = 28, y = 510, size = 12.5)
)

run_figure(cfg, edits)

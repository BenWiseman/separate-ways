#!/usr/bin/env Rscript

script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1]])
source(file.path(dirname(normalizePath(script_file)), "_helpers.R"))

cfg <- list(
  name = "fig1_dictionary",
  pdf_width_pt = 504,
  pdf_height_pt = 417.6,
  png_window_width = 1400,
  png_window_height = 1160,
  png_scale = 1,
  png_use_wrapper = TRUE
)

# These are the labels most likely to need adjustment after copy edits. Run
# with --list-text to print every text node and its current coordinates.
edits <- list(
  text_edit("Figure 1. Where a fold can sit", x = 28, y = 29, size = 17),
  text_edit("Shaded rows are consistent", x = 28, y = 48, size = 12.5),
  text_edit("fixed set", x = 30, y = 78, size = 13, exact = TRUE),
  text_edit("identification", x = 138, y = 78, size = 13, exact = TRUE),
  text_edit("price", x = 321, y = 78, size = 13, exact = TRUE),
  text_edit("verdict", x = 551, y = 78, size = 13, exact = TRUE),
  text_edit("In the two-copy row", x = 28, y = 563, size = 12.5),

  # --- legibility fixes 2026-09-17 ---
  # The longest bold identification is wrapped at the full 12.5-unit body
  # size. Shrinking it would fall below 8 pt after the manuscript scales the
  # seven-inch source to its 6.5-inch text block.
  text_edit("glue at the bifurcation sphere", x = 141, y = 459, size = 12.5,
            lines = c("glue at the", "bifurcation sphere"), line_height = 15,
            exact = TRUE),
  text_edit("second sheet is the J-image", y = 497, exact = TRUE),

  # Column 1 is centred on x = 78 but its
  # header hung left-anchored at x = 30, off the column's optical axis.
  text_edit("fixed set", x = 78, anchor = "middle", exact = TRUE),

  # Row captions are column-1 sub-labels, not signals. Bold made them compete
  # with the identification column. Medium weight, same size, darker fill.
  text_edit("no crease", weight = 500, fill = "#1a1a1a", exact = TRUE),
  text_edit("codimension 1", weight = 500, fill = "#1a1a1a", exact = TRUE, occurrence = 1),
  text_edit("codimension 2", weight = 500, fill = "#1a1a1a", exact = TRUE),
  text_edit("codimension 1", weight = 500, fill = "#1a1a1a", exact = TRUE, occurrence = 2),
  text_edit("two copies", weight = 500, fill = "#1a1a1a", exact = TRUE),

  # Icons now end at T+56, so captions drop to the baseline T+74.
  # Gap icon-to-cap-top is 9px, clearance to the row rule below is 11px.
  text_edit("no crease", y = 162, exact = TRUE),
  text_edit("codimension 1", y = 250, exact = TRUE, occurrence = 1),
  text_edit("codimension 2", y = 338, exact = TRUE),
  text_edit("codimension 1", y = 426, exact = TRUE, occurrence = 2),
  text_edit("two copies", y = 522, exact = TRUE),

  # Grey secondary text over hatching was the real legibility failure.
  text_edit("fixes the bifurcation 2-sphere", fill = "#333333", exact = TRUE),
  text_edit("+\u2212+ three-brane model", fill = "#333333", exact = TRUE)

)

run_figure(cfg, edits)

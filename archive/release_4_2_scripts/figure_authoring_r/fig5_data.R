#!/usr/bin/env Rscript

script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1]])
source(file.path(dirname(normalizePath(script_file)), "_helpers.R"))

cfg <- list(
  name = "fig5_data",
  pdf_width_pt = 504,
  pdf_height_pt = 630,
  png_window_width = 672,
  png_window_height = 840,
  png_scale = 2,
  png_use_wrapper = FALSE
)

edits <- list(
  text_edit("(a) Neutrino-mass markers", x = 18, y = 18, size = "11pt"),
  text_edit("circles: exact-rule markers", x = 18, y = 33, size = "10pt"),
  text_edit("(b) Geometry-only mass profiles", x = 18, y = 260, size = "11pt"),
  text_edit("BAO + Pantheon+", x = 18, y = 274, size = "10pt"),
  text_edit("Each profile", x = 18, y = 288, size = "10pt"),
  text_edit("(c) Exact-Λ stance", x = 18, y = 470, size = "11pt"),
  # Priority 0 from the 180dpi audit, confirmed by measuring: this label
  # sat at x = 135.67 with anchor "end", so it ran LEFTWARD across the 0.2 y
  # tick at x = 67 and into the top plot border. Re-anchored to start and
  # moved below-right of the (-1, 0) marker, into clear space. The marker
  # itself is untouched and still sits at exactly (-1, 0).
  text_edit("exact \u039b: (\u22121, 0)", x = 155, y = 521,
            anchor = "start", size = "10pt"),
  text_edit("3 Unite; bars", x = 310, y = 526, size = "10pt"),

  # Panel (a) now uses two semantic groups and five single-line rows. Retain
  # 10 pt labels so final manuscript type is 8.03 pt at 5.6225 inches wide.
  text_edit("EXACT RULE", x = 18, y = 59, size = "10pt", exact = TRUE),
  text_edit("NO · m\u2081 = 0", x = 18, y = 82, size = "10pt", exact = TRUE),
  text_edit("IO · m\u2083 = 0", x = 18, y = 105, size = "10pt", exact = TRUE),
  text_edit("PUBLISHED LIMITS", x = 18, y = 134, size = "10pt", exact = TRUE),
  text_edit("Feldman\u2013Cousins", size = "10pt", exact = TRUE),
  text_edit("DESI DR2 · \u039bCDM", size = "10pt", exact = TRUE),
  text_edit("DESI DR2 · w\u2080w\u2090CDM", size = "10pt", exact = TRUE)
)

run_figure(cfg, edits)

#!/usr/bin/env Rscript

script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1]])
source(file.path(dirname(normalizePath(script_file)), "_helpers.R"))

cfg <- list(
  name = "graphical_abstract_v3",
  pdf_width_pt = 1200,
  pdf_height_pt = 855,
  png_window_width = 1600,
  png_window_height = 1140,
  png_scale = 2,
  png_use_wrapper = FALSE
)

# Adjustment panel for the spacious two-zone composition. The cosmological
# drawing is explicitly interpretive; the lower card contains the demonstrated
# causally disjoint free-field construction. Geometry stays in the SVG master.
edits <- list(
  text_edit("Two CPT-related branches, one retained quantum algebra",
            x = 70, y = 62, size = 42, exact = TRUE),
  text_edit("INTERPRETIVE SCHEMATIC", x = 1327.5, y = 105, size = 28,
            anchor = "middle", weight = 700, exact = TRUE),
  text_edit("CPT-related branch", x = 345, y = 150, size = 31, exact = TRUE),
  text_edit("CPT-symmetric boundary", x = 800, y = 142, size = 31, exact = TRUE),
  text_edit("our branch", x = 1255, y = 150, size = 31, exact = TRUE),
  text_edit("local future", x = 455, y = 188, size = 29, occurrence = 1),
  text_edit("local future", x = 1145, y = 188, size = 29, occurrence = 2),
  text_edit("Θ = CPT", x = 800, y = 272, size = 29, exact = TRUE),
  text_edit("fixed reference future", x = 800, y = 315, size = 28),
  text_edit("time-reversed portrait", x = 360, y = 505, size = 28),
  text_edit("white-hole orientation (global)", x = 360, y = 538, size = 28),
  text_edit("future-absorbing horizon", x = 1240, y = 505, size = 28),
  text_edit("black-hole orientation", x = 1240, y = 538, size = 28),
  text_edit("charge-conjugate, parity-reflected", x = 360, y = 662, size = 29),
  text_edit("ordinary matter orientation", x = 1240, y = 662, size = 29),

  text_edit("DERIVED", x = 235, y = 771, size = 28, fill = "#ffffff",
            weight = 700, exact = TRUE),
  text_edit("Shared average/difference algebra", x = 345, y = 771, size = 34),
  text_edit("stated causally disjoint free-field pairing", x = 345, y = 806, size = 28),
  text_edit("average field", x = 460, y = 862, size = 30, exact = TRUE),
  text_edit("Φc = ½(Φ + ΘΦ)", x = 460, y = 908, size = 36,
            subscript_size = 30, subscript_shift = "-22%"),
  text_edit("classical Gaussian marginal", x = 460, y = 948, size = 29),
  text_edit("difference field", x = 1140, y = 862, size = 30, exact = TRUE),
  text_edit("Φq = Φ − ΘΦ", x = 1140, y = 908, size = 36,
            subscript_size = 30, subscript_shift = "-22%"),
  text_edit("retains the quantum conjugate", x = 1140, y = 948, size = 29),
  text_edit("[Φc, Φc]", x = 625, y = 1018, size = 31,
            subscript_size = 28, subscript_shift = "-22%"),
  text_edit("OPEN · interactions", x = 1287.5, y = 1017, size = 28,
            anchor = "middle", weight = 700, exact = TRUE),
  text_edit("The drawing is interpretive", x = 800, y = 1120, size = 28)
)

run_figure(cfg, edits)

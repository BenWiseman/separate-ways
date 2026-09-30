#!/usr/bin/env Rscript

script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1]])
source(file.path(dirname(normalizePath(script_file)), "_helpers.R"))

cfg <- list(
  name = "fig2_janus",
  pdf_width_pt = 504,
  pdf_height_pt = 360,
  png_window_width = 1400,
  png_window_height = 1000,
  png_scale = 1,
  png_use_wrapper = TRUE
)

edits <- list(
  text_edit("Figure 2. The bang as a complexity minimum", x = 28, y = 29, size = 17),
  text_edit("Positive Bel", x = 28, y = 48, size = 12.5),
  text_edit("Regular CPT modes", x = 28, y = 72, size = 12.5),
  # Keep the panel tag clear of the two descending irregular-mode curves.
  text_edit("(a) scaling", x = 312, y = 154, size = 13, anchor = "end"),
  text_edit("(b) regular modes", x = 402, y = 154, size = 13),
  text_edit("bang, η = 0", x = 518, y = 174, size = 12.5),
  text_edit("anti-verse", x = 403, y = 336, size = 12.5),
  text_edit("this-verse", x = 633, y = 336, size = 12.5),
  text_edit("conformal time η", x = 518, y = 394, size = 13, occurrence = 2),
  text_edit("P/Pmax is exactly even", x = 518, y = 416, size = 12.5),
  text_edit("Regular modes vanish", x = 28, y = 444, size = 12.5),
  text_edit("Linear order only", x = 28, y = 466, size = 12.5),
  # The "bang, eta = 0" label sat at y = 174 with the dashed vertical line
  # starting at y = 175, both centred on x = 518, so the line emerged through
  # the middle of the words. Label lifted to 168, line start moved to 184 in
  # the SVG, giving a clear 16px gap.
  text_edit("bang, \u03b7 = 0", y = 168, exact = TRUE)
)

# The exact points in the release SVG come from compute_figs.py.  The original
# compact reflow used a vertical range whose top coincided with the irregular
# tensor curve and a horizontal range ending at the final x tick.  Re-map those
# same points into a panel with a half-decade of right-side space and a 10-unit
# top margin.  A clip path remains a guard against later line-width changes.
remap_left_panel <- function(svg) {
  doc <- xml2::read_xml(svg)
  curves <- xml2::xml_find_all(doc, ".//*[local-name()='polyline']")
  if (length(curves) < 4L) stop("Figure 2 is missing its four left-panel curves.")
  if (!is.na(xml2::xml_attr(curves[[1]], "clip-path"))) return(invisible(FALSE))

  remap_points <- function(points) {
    pairs <- strsplit(trimws(points), "[[:space:]]+")[[1]]
    mapped <- vapply(pairs, function(pair) {
      xy <- as.numeric(strsplit(pair, ",", fixed = TRUE)[[1]])
      # x: log10|eta| [-5, -2] -> [-5, -1.4]; y: [-9, 2.5] -> [-9, 3.5].
      x <- 58 + (xy[[1]] - 58) * (3 / 3.6)
      y <- 28 + 0.92 * xy[[2]]
      sprintf("%.2f,%.2f", x, y)
    }, character(1))
    paste(mapped, collapse = " ")
  }
  for (i in 1:4) {
    xml2::xml_set_attr(curves[[i]], "points", remap_points(xml2::xml_attr(curves[[i]], "points")))
  }

  # Move the left-panel y ticks to the expanded logarithmic range.
  y_axis_ticks <- xml2::xml_find_all(doc,
    ".//*[local-name()='line' and @x1='53.0' and @x2='58.0']")
  for (tick in y_axis_ticks) {
    for (attr in c("y1", "y2")) {
      old <- as.numeric(xml2::xml_attr(tick, attr))
      xml2::xml_set_attr(tick, attr, sprintf("%.2f", 28 + 0.92 * old))
    }
  }
  y_axis_labels <- xml2::xml_find_all(doc,
    ".//*[local-name()='text' and @x='48.0']")
  for (label in y_axis_labels) {
    old <- as.numeric(xml2::xml_attr(label, "y"))
    xml2::xml_set_attr(label, "y", sprintf("%.2f", 28 + 0.92 * old))
  }

  # The x ticks remain the same values but gain the restored right-hand range.
  x_ticks <- xml2::xml_find_all(doc,
    ".//*[local-name()='line' and @y1='350.0' and @y2='355.0' and @x1=@x2]")
  for (tick in x_ticks) {
    old <- as.numeric(xml2::xml_attr(tick, "x1"))
    if (old <= 328) {
      new <- 58 + (old - 58) * (3 / 3.6)
      xml2::xml_set_attr(tick, "x1", sprintf("%.2f", new))
      xml2::xml_set_attr(tick, "x2", sprintf("%.2f", new))
    }
  }
  x_tick_labels <- xml2::xml_find_all(doc,
    ".//*[local-name()='text' and @y='370.0']")
  for (label in x_tick_labels) {
    old <- as.numeric(xml2::xml_attr(label, "x"))
    if (old <= 328) xml2::xml_set_attr(label, "x", sprintf("%.2f", 58 + (old - 58) * (3 / 3.6)))
  }

  defs <- xml2::xml_find_first(doc, ".//*[local-name()='defs']")
  clip <- xml2::xml_add_child(defs, "clipPath")
  xml2::xml_set_attr(clip, "id", "fig2-left-plot-clip")
  clip_rect <- xml2::xml_add_child(clip, "rect")
  xml2::xml_set_attrs(clip_rect, c(x = "58", y = "92", width = "270", height = "258"))
  for (i in 1:4) xml2::xml_set_attr(curves[[i]], "clip-path", "url(#fig2-left-plot-clip)")

  xml2::write_xml(doc, svg)
  svg_text <- paste(readLines(svg, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  svg_text <- gsub("</tspan>[[:space:]]+<tspan", "</tspan><tspan", svg_text, perl = TRUE)
  writeLines(svg_text, svg, useBytes = TRUE)
  invisible(TRUE)
}

cli <- parse_cli()
result <- run_figure(cfg, edits)
if (remap_left_panel(result$svg) && !cli$no_render) {
  render_svg(result$svg, result$pdf, result$png, cfg)
}

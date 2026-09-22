# Shared authoring helpers for the Paper 2 V3 SVG figures.
#
# The release SVG remains the exact-content base. Each figure script applies
# named edits to a copy and writes a preview under tmp/ unless --publish is
# passed. This keeps exploratory layout work away from submission assets.

if (!requireNamespace("xml2", quietly = TRUE)) {
  stop("The CRAN package 'xml2' is required. Run install.packages('xml2').")
}

authoring_dir <- local({
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (!length(file_arg)) stop("Run these files with Rscript.")
  dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = TRUE))
})

project_root <- normalizePath(file.path(authoring_dir, "../../.."), mustWork = TRUE)
live_dir <- file.path(project_root, "pub", "paper2")
preview_root <- file.path(project_root, "tmp", "figure_authoring_r")

normalise_text <- function(x) {
  trimws(gsub("[[:space:]]+", " ", x))
}

text_edit <- function(match, x = NULL, y = NULL, size = NULL,
                      anchor = NULL, weight = NULL, fill = NULL,
                      value = NULL, occurrence = 1L, exact = FALSE,
                      max_width = NULL, subscript_size = NULL,
                      subscript_shift = NULL, lines = NULL,
                      line_height = NULL) {
  list(
    match = match, x = x, y = y, size = size, anchor = anchor,
    weight = weight, fill = fill, value = value,
    occurrence = as.integer(occurrence), exact = exact,
    max_width = max_width, subscript_size = subscript_size,
    subscript_shift = subscript_shift, lines = lines,
    line_height = line_height
  )
}

parse_cli <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  out_arg <- grep("^--out-dir=", args, value = TRUE)
  list(
    publish = "--publish" %in% args,
    list_text = "--list-text" %in% args,
    no_render = "--no-render" %in% args,
    out_dir = if (length(out_arg)) sub("^--out-dir=", "", out_arg[[1]]) else NULL
  )
}

svg_text_nodes <- function(doc) {
  xml2::xml_find_all(doc, ".//*[local-name()='text']")
}

list_text_nodes <- function(doc) {
  nodes <- svg_text_nodes(doc)
  rows <- lapply(seq_along(nodes), function(i) {
    node <- nodes[[i]]
    data.frame(
      index = i,
      text = normalise_text(xml2::xml_text(node)),
      x = xml2::xml_attr(node, "x") %||% "",
      y = xml2::xml_attr(node, "y") %||% "",
      size = xml2::xml_attr(node, "font-size") %||% "",
      class = xml2::xml_attr(node, "class") %||% "",
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  print(out, row.names = FALSE, right = FALSE)
  invisible(out)
}

`%||%` <- function(x, y) if (is.na(x) || !length(x)) y else x

find_text_node <- function(doc, edit) {
  nodes <- svg_text_nodes(doc)
  labels <- vapply(nodes, function(node) normalise_text(xml2::xml_text(node)), character(1))
  hits <- if (isTRUE(edit$exact)) {
    which(labels == normalise_text(edit$match))
  } else {
    which(grepl(normalise_text(edit$match), labels, fixed = TRUE))
  }
  if (!length(hits)) {
    stop("No SVG text node matched: ", edit$match)
  }
  if (edit$occurrence < 1L || edit$occurrence > length(hits)) {
    stop("Occurrence ", edit$occurrence, " is out of range for: ", edit$match)
  }
  nodes[[hits[[edit$occurrence]]]]
}

set_optional_attr <- function(node, name, value) {
  if (!is.null(value)) xml2::xml_set_attr(node, name, as.character(value))
}

apply_text_edit <- function(doc, edit) {
  node <- find_text_node(doc, edit)

  if (!is.null(edit$value) && !is.null(edit$lines)) {
    stop("Use either value or lines in a text edit, not both: ", edit$match)
  }

  # Setting value replaces child tspans. Omit value when the equation or label
  # already has structured tspans that should remain intact.
  if (!is.null(edit$value)) xml2::xml_set_text(node, edit$value)

  # Wrap a long label without reducing its type size. A trailing literal
  # space on every non-final line preserves the original normalised text, so
  # the edit remains idempotent on later authoring passes.
  if (!is.null(edit$lines)) {
    if (length(edit$lines) < 2L || is.null(edit$line_height)) {
      stop("Wrapped text needs at least two lines and line_height: ", edit$match)
    }
    xml2::xml_remove(xml2::xml_children(node))
    xml2::xml_set_text(node, "")
    line_x <- edit$x %||% xml2::xml_attr(node, "x")
    for (i in seq_along(edit$lines)) {
      value <- paste0(edit$lines[[i]], if (i < length(edit$lines)) " " else "")
      span <- xml2::xml_add_child(node, "tspan", value)
      xml2::xml_set_attr(span, "x", as.character(line_x))
      xml2::xml_set_attr(span, "dy", if (i == 1L) "0" else as.character(edit$line_height))
    }
  }

  set_optional_attr(node, "x", edit$x)
  set_optional_attr(node, "y", edit$y)
  set_optional_attr(node, "font-size", edit$size)
  set_optional_attr(node, "text-anchor", edit$anchor)
  set_optional_attr(node, "font-weight", edit$weight)
  set_optional_attr(node, "fill", edit$fill)

  # max_width is an explicit last-resort fit control. It preserves the chosen
  # font size and lets SVG adjust spacing. Moving or wrapping text is usually
  # more readable, so no script enables this by default.
  if (!is.null(edit$max_width)) {
    xml2::xml_set_attr(node, "textLength", as.character(edit$max_width))
    xml2::xml_set_attr(node, "lengthAdjust", "spacingAndGlyphs")
  }

  if (!is.null(edit$subscript_size) || !is.null(edit$subscript_shift)) {
    # Match both untouched SVG subscripts ("sub") and an explicit percentage
    # written by an earlier authoring pass, so the edit is idempotent.
    tspans <- xml2::xml_find_all(node, ".//*[local-name()='tspan' and @baseline-shift]")
    if (!length(tspans)) stop("No subscript tspans found in: ", edit$match)
    if (!is.null(edit$subscript_size)) {
      xml2::xml_set_attr(tspans, "font-size", as.character(edit$subscript_size))
    }
    if (!is.null(edit$subscript_shift)) {
      xml2::xml_set_attr(tspans, "baseline-shift", as.character(edit$subscript_shift))
    }
  }

  invisible(node)
}

apply_text_edits <- function(doc, edits) {
  if (!length(edits)) return(invisible(doc))
  for (edit in edits) apply_text_edit(doc, edit)
  invisible(doc)
}

chrome_run <- function(args) {
  chrome <- Sys.which("google-chrome")
  if (!nzchar(chrome)) stop("google-chrome is required for PDF and PNG previews.")
  result <- suppressWarnings(system2(chrome, c(paste0("--user-data-dir=", file.path(tempdir(), "chrome-figure")), args), stdout = TRUE, stderr = TRUE))
  status <- attr(result, "status") %||% 0L
  if (status != 0L) stop(paste(result, collapse = "\n"))
  invisible(result)
}

render_svg <- function(svg, pdf, png, cfg) {
  render_dir <- file.path(preview_root, "render_html")
  dir.create(render_dir, recursive = TRUE, showWarnings = FALSE)
  html <- tempfile(pattern = paste0(cfg$name, "_pdf_"), tmpdir = render_dir, fileext = ".html")
  png_html <- tempfile(pattern = paste0(cfg$name, "_png_"), tmpdir = render_dir, fileext = ".html")
  svg_uri <- paste0("file://", normalizePath(svg, mustWork = TRUE))
  html_text <- sprintf(
    paste0(
      '<style>@page{size:%spt %spt;margin:0}',
      'html,body,img{margin:0;width:%spt;height:%spt;display:block;overflow:hidden}</style>',
      '<img src="%s">'
    ),
    cfg$pdf_width_pt, cfg$pdf_height_pt,
    cfg$pdf_width_pt, cfg$pdf_height_pt, svg_uri
  )
  writeLines(html_text, html, useBytes = TRUE)
  png_html_text <- sprintf(
    paste0(
      '<style>html,body,img{margin:0;width:%spx;height:%spx;',
      'display:block;overflow:hidden}</style><img src="%s">'
    ),
    cfg$png_window_width, cfg$png_window_height, svg_uri
  )
  writeLines(png_html_text, png_html, useBytes = TRUE)

  chrome_run(c(
    "--headless", "--no-sandbox", "--disable-gpu", "--no-pdf-header-footer",
    paste0("--print-to-pdf=", normalizePath(dirname(pdf), mustWork = TRUE), "/", basename(pdf)),
    paste0("file://", html)
  ))
  chrome_run(c(
    "--headless", "--no-sandbox", "--disable-gpu", "--hide-scrollbars",
    paste0("--force-device-scale-factor=", cfg$png_scale),
    paste0("--window-size=", cfg$png_window_width, ",", cfg$png_window_height),
    paste0("--screenshot=", normalizePath(dirname(png), mustWork = TRUE), "/", basename(png)),
    if (isTRUE(cfg$png_use_wrapper)) paste0("file://", png_html) else svg_uri
  ))
  unlink(c(html, png_html))
}

backup_release_assets <- function(name) {
  stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  backup_dir <- file.path(preview_root, "backups", paste0(name, "_", stamp))
  dir.create(backup_dir, recursive = TRUE, showWarnings = FALSE)
  files <- file.path(live_dir, paste0(name, c(".svg", ".pdf", ".png")))
  files <- files[file.exists(files)]
  if (length(files) && !all(file.copy(files, backup_dir, overwrite = FALSE))) {
    stop("Could not back up the current release assets.")
  }
  backup_dir
}

run_figure <- function(cfg, edits = list()) {
  cli <- parse_cli()
  input_svg <- file.path(live_dir, paste0(cfg$name, ".svg"))
  if (!file.exists(input_svg)) stop("Missing base SVG: ", input_svg)
  doc <- xml2::read_xml(input_svg)

  if (cli$list_text) {
    list_text_nodes(doc)
    return(invisible(NULL))
  }

  apply_text_edits(doc, edits)

  out_dir <- if (cli$publish) {
    backup <- backup_release_assets(cfg$name)
    message("Backed up release assets to: ", backup)
    live_dir
  } else if (!is.null(cli$out_dir)) {
    normalizePath(cli$out_dir, mustWork = FALSE)
  } else {
    file.path(preview_root, cfg$name)
  }
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  svg <- file.path(out_dir, paste0(cfg$name, ".svg"))
  pdf <- file.path(out_dir, paste0(cfg$name, ".pdf"))
  png <- file.path(out_dir, paste0(cfg$name, ".png"))
  # xml2 formats child elements on separate lines. In SVG text, whitespace
  # between adjacent tspans becomes visible and breaks equation spacing. Join
  # those boundaries again after serialization.
  xml2::write_xml(doc, svg)
  svg_text <- paste(readLines(svg, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  svg_text <- gsub("</tspan>[[:space:]]+<tspan", "</tspan><tspan", svg_text, perl = TRUE)
  writeLines(svg_text, svg, useBytes = TRUE)

  if (!cli$no_render) render_svg(svg, pdf, png, cfg)

  message("SVG: ", svg)
  if (!cli$no_render) {
    message("PDF: ", pdf)
    message("PNG: ", png)
  }
  if (cli$publish) {
    message("Release assets changed. Rebuild the manuscript PDFs and SHA256SUMS_v3.txt before submission.")
  }
  invisible(list(svg = svg, pdf = pdf, png = png))
}

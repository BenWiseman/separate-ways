#!/usr/bin/env Rscript
# Does any label a figure draws land on top of something else the figure drew?
#
#   Rscript checks/label_ink_check.R checks/fig_*.R
#   Rscript checks/label_ink_check.R --selftest
#
# WHY. label_fit_check.R measures a label against its PANEL and catches text that runs off the
# edge. It says nothing about what is already inside the panel. On 2026-09-28 the companion was
# read as a PDF and Figure 5's "double-precision floor" was printed over by the bars at l = 5, 7
# and 9: inside the panel, correctly sized, and unreadable. Gate 19 passed it. This is the check
# that would not have.
#
# HOW. Source each generator with text() and the drawing primitives shadowed, collecting label
# boxes and ink in user coordinates, then test every label box against every piece of ink drawn
# in the same panel. Panels are told apart by par("fig"), so a label in one panel is never
# compared with ink in another.
#
# WHAT IT DOES NOT COVER, said plainly so the pass means something: axis(), box() and the tick
# labels R draws itself are not shadowed, because a label sitting on an axis line is not a
# defect. A rect filling more than half its panel is treated as a background and skipped, for
# the same reason. Text boxes are shrunk to 70 per cent before testing, so a label that merely
# abuts a line is spared and only one genuinely written through is reported.

.probe.args <- commandArgs(TRUE)
.ink <- new.env()
.ink$labels <- list(); .ink$marks <- list(); .ink$covers <- list()
.ink$hits <- character(0); .ink$ok <- character(0); .ink$name <- "?"

.panel <- function() paste(signif(par("fig"), 6), collapse = ",")
# THE TRAP label_fit_check.R's header already warns about, walked into anyway on the first run:
# strwidth and strheight come back in USER coordinates, which are log10 units on a log axis,
# while the x and y handed to text() and lines() are raw data. Mixing them put a label at y = 0
# and the line that crosses it at y = 3e-4, and the checker reported nonsense collisions on
# every log-axis figure. Everything is converted on the way in now.
# Ink that cannot hide text is not ink. A pale shaded band or a pale guide line leaves every
# letter legible, and flagging them buries the two or three labels that are genuinely lost. The
# test is the fill's own relative luminance, so it needs no list of colours to maintain.
.pale <- function(col) {
  if (is.null(col)) return(FALSE)          # no colour given means the default, which is black
  if (all(is.na(col))) return(TRUE)
  v <- tryCatch(grDevices::col2rgb(col[1])/255, error = function(e) matrix(0, 3, 1))
  as.numeric(0.2126*v[1] + 0.7152*v[2] + 0.0722*v[3]) > 0.80
}
.ux <- function(x) if (par("xlog")) suppressWarnings(log10(x)) else x
.uy <- function(y) if (par("ylog")) suppressWarnings(log10(y)) else y
.usr2  <- function() { u <- par("usr"); c(min(u[1],u[2]), max(u[1],u[2]), min(u[3],u[4]), max(u[3],u[4])) }

text <- function(x, y = NULL, labels = seq_along(x), ...) {
  o <- list(...)
  cex <- if (!is.null(o$cex)) o$cex else 1
  adj <- if (!is.null(o$adj)) o$adj[1] else 0.5
  lab <- if (is.expression(labels) || is.call(labels)) character(0) else as.character(labels)
  if (length(lab) && is.null(o$srt)) {
    for (i in seq_along(lab)) {
      xi <- .ux(x[[min(i, length(x))]]); yi <- if (is.null(y)) NA else .uy(y[[min(i, length(y))]])
      if (is.finite(xi) && is.finite(yi)) {
        w <- abs(strwidth(lab[i], cex = cex)); h <- abs(strheight(lab[i], cex = cex))
        x0 <- xi - adj*w; x1 <- x0 + w; y0 <- yi - 0.5*h; y1 <- yi + 0.5*h
        .ink$labels[[length(.ink$labels) + 1]] <-
          list(p = .panel(), lab = lab[i], b = c(min(x0,x1), max(x0,x1), y0, y1), u = .usr2())
      }
    }
  }
  graphics::text(x, y, labels, ...)
}
.addseg <- function(xs, ys, col = NULL) {
  if (.pale(col)) return(invisible())
  xs <- .ux(as.numeric(xs)); ys <- .uy(as.numeric(ys))
  k <- min(length(xs), length(ys)); if (k < 2) return(invisible())
  ok <- is.finite(xs[1:k]) & is.finite(ys[1:k])
  for (i in seq_len(k - 1)) if (ok[i] && ok[i+1])
    .ink$marks[[length(.ink$marks) + 1]] <-
      list(p = .panel(), kind = "seg", v = c(xs[i], ys[i], xs[i+1], ys[i+1]))
  invisible()
}
lines    <- function(x, y = NULL, ...) { if (is.null(y) && is.list(x)) { y <- x$y; x <- x$x }
                                         .addseg(x, y, list(...)$col); graphics::lines(x, y, ...) }
polygon  <- function(x, y = NULL, ...) { .addseg(c(x, x[1]), c(y, y[1]), list(...)$col); graphics::polygon(x, y, ...) }
points   <- function(x, y = NULL, ...) { .addseg(c(x, x), c(y, y), list(...)$col); graphics::points(x, y, ...) }
segments <- function(x0, y0, x1 = x0, y1 = y0, ...) {
  n <- max(length(x0), length(y0), length(x1), length(y1))
  for (i in seq_len(n)) .addseg(c(x0[[min(i,length(x0))]], x1[[min(i,length(x1))]]),
                                c(y0[[min(i,length(y0))]], y1[[min(i,length(y1))]]), list(...)$col)
  graphics::segments(x0, y0, x1, y1, ...)
}
arrows   <- function(x0, y0, x1 = x0, y1 = y0, ...) {
  n <- max(length(x0), length(y0), length(x1), length(y1))
  for (i in seq_len(n)) .addseg(c(x0[[min(i,length(x0))]], x1[[min(i,length(x1))]]),
                                c(y0[[min(i,length(y0))]], y1[[min(i,length(y1))]]), list(...)$col)
  graphics::arrows(x0, y0, x1, y1, ...)
}
abline   <- function(a = NULL, b = NULL, h = NULL, v = NULL, ...) {
  u <- .usr2()
  if (!.pale(list(...)$col)) {
    for (hh in h) { yy <- .uy(hh); .ink$marks[[length(.ink$marks)+1]] <-
                      list(p = .panel(), kind = "seg", v = c(u[1], yy, u[2], yy)) }
    for (vv in v) { xx <- .ux(vv); .ink$marks[[length(.ink$marks)+1]] <-
                      list(p = .panel(), kind = "seg", v = c(xx, u[3], xx, u[4])) }
  }
  graphics::abline(a, b, h, v, ...)
}
rect <- function(xleft, ybottom, xright, ytop, ...) {
  o <- list(...); filled <- !is.null(o$col) && !all(is.na(o$col)) && !.pale(o$col)
  u <- .usr2(); area <- (u[2]-u[1])*(u[4]-u[3])
  n <- max(length(xleft), length(ybottom), length(xright), length(ytop))
  for (i in seq_len(n)) {
    xl <- .ux(xleft[[min(i,length(xleft))]]); xr <- .ux(xright[[min(i,length(xright))]])
    yb <- .uy(ybottom[[min(i,length(ybottom))]]); yt <- .uy(ytop[[min(i,length(ytop))]])
    if (!all(is.finite(c(xl,xr,yb,yt)))) next
    frac <- abs((xr-xl)*(yt-yb))/max(area, .Machine$double.eps)
    # A shaded HALF of a panel is a background too, and a band spanning a whole axis is one
    # whatever its area, so neither counts as ink a label could be lost in.
    spanx <- abs(xr-xl) >= 0.95*(u[2]-u[1]); spany <- abs(yt-yb) >= 0.95*(u[4]-u[3])
    box <- c(min(xl,xr), max(xl,xr), min(yb,yt), max(yb,yt))
    # A small filled rect can MASK what is behind it, pale or not, and that is how a label is
    # rescued where it has to stay put. A panel-sized one is the background and masks nothing,
    # so it must not count, or every label in the figure is spared.
    if (!is.null(o$col) && !all(is.na(o$col)) && frac < 0.45 && !spanx && !spany)
      .ink$covers[[length(.ink$covers) + 1]] <- list(p = .panel(), v = box)
    if (filled && frac < 0.45 && !spanx && !spany)
      .ink$marks[[length(.ink$marks) + 1]] <- list(p = .panel(), kind = "box", v = box)
    else if (is.null(o$col) || all(is.na(o$col))) .addseg(c(xl,xr,xr,xl,xl), c(yb,yb,yt,yt,yb), o$border)
  }
  graphics::rect(xleft, ybottom, xright, ytop, ...)
}

# segment against axis-aligned box, Liang-Barsky
.seg_box <- function(s, b) {
  dx <- s[3]-s[1]; dy <- s[4]-s[2]; t0 <- 0; t1 <- 1
  for (k in 1:4) {
    p <- c(-dx, dx, -dy, dy)[k]
    q <- c(s[1]-b[1], b[2]-s[1], s[2]-b[3], b[4]-s[2])[k]
    if (p == 0) { if (q < 0) return(FALSE) } else {
      r <- q/p
      if (p < 0) { if (r > t1) return(FALSE); if (r > t0) t0 <- r }
      else       { if (r < t0) return(FALSE); if (r < t1) t1 <- r }
    }
  }
  TRUE
}
.box_box <- function(a, b) a[1] < b[2] && a[2] > b[1] && a[3] < b[4] && a[4] > b[3]

.report <- function() {
  for (L in .ink$labels) {
    b <- L$b
    cx <- 0.5*(b[1]+b[2]); cy <- 0.5*(b[3]+b[4])
    s <- 0.70                      # shrink, so abutting is spared and only crossing is reported
    bb <- c(cx - s*(cx-b[1]), cx + s*(b[2]-cx), cy - s*(cy-b[3]), cy + s*(b[4]-cy))
    # A filled shape that CONTAINS the label is not a defect: it is either a labelled box, which
    # is how every flow chart in this release is drawn, or a patch put there on purpose to mask
    # what is behind the text. Either way the label reads. Only a shape that covers PART of it,
    # or a line through it, is the thing this check is for.
    covered <- FALSE
    for (M in .ink$covers)
      if (M$p == L$p &&
          M$v[1] <= bb[1] && M$v[2] >= bb[2] && M$v[3] <= bb[3] && M$v[4] >= bb[4]) {
        covered <- TRUE; break
      }
    if (covered) next
    for (M in .ink$marks) {
      if (M$p != L$p) next
      hit <- if (M$kind == "seg") .seg_box(M$v, bb) else .box_box(M$v, bb)
      if (hit) {
        where <- if (M$kind == "seg")
          sprintf("a line from (%.4g, %.4g) to (%.4g, %.4g)", M$v[1], M$v[2], M$v[3], M$v[4])
        else sprintf("a filled shape spanning x %.4g..%.4g, y %.4g..%.4g", M$v[1], M$v[2], M$v[3], M$v[4])
        if (!.is.accepted(.ink$name, L$lab))
          .ink$hits <- c(.ink$hits, sprintf('%s: "%s" at (%.4g, %.4g) is drawn over by %s',
                                            .ink$name, L$lab, 0.5*(b[1]+b[2]), 0.5*(b[3]+b[4]), where))
        else .ink$ok <- c(.ink$ok, sprintf("%s: \"%s\"", .ink$name, L$lab))
        break
      }
    }
  }
}

# Accepted overlaps, with the reason, in the file beside this one. A label on that list is
# reported as accepted rather than as a defect, so a NEW overlap still fails while the four
# that were looked at and left alone do not.
.accept <- list()
.acc.f <- "checks/ink_overlaps_accepted.tsv"
if (file.exists(.acc.f)) {
  for (.ln in readLines(.acc.f)) {
    if (!nzchar(trimws(.ln)) || substr(trimws(.ln), 1, 1) == "#") next
    .p <- strsplit(.ln, "\t")[[1]]
    if (length(.p) >= 2) .accept[[length(.accept) + 1]] <- c(.p[1], .p[2])
  }
}
.is.accepted <- function(gen, lab) {
  for (a in .accept) if (a[1] == gen && a[2] == lab) return(TRUE)
  FALSE
}

if ("--selftest" %in% .probe.args) {
  pdf(tempfile(fileext = ".pdf"), width = 5, height = 4)
  .ink$name <- "plant"
  plot(NA, xlim = c(0, 10), ylim = c(0, 10), axes = FALSE, xlab = "", ylab = "")
  rect(0, 0, 10, 10, col = "#eeeeee")         # the panel background, drawn first as figures do
  segments(0, 5, 10, 5, col = "black")        # a line straight through the label below
  text(5, 5, "written over")                  # must be caught
  text(5, 9, "clear of everything")           # must be spared
  text(5, 2, "on the background only")        # must be spared
  rect(3.4, 7.4, 6.6, 7.9, col = "#eeeeee")   # a mask under the label below
  segments(0, 7.65, 10, 7.65, col = "black")  # a line the mask protects it from
  text(5, 7.65, "masked")                     # must be spared
  dev.off()
  .report()
  caught <- any(grepl("written over", .ink$hits))
  spared3 <- !any(grepl("masked", .ink$hits))
  spared1 <- !any(grepl("clear of everything", .ink$hits))
  spared2 <- !any(grepl("on the background only", .ink$hits))
  cat(sprintf("   plant: a label crossed by a line is caught: %s\n", if (caught) "yes" else "NO"))
  cat(sprintf("   plant: a label clear of everything is spared: %s\n", if (spared1) "yes" else "NO"))
  cat(sprintf("   plant: a label on a panel background is spared: %s\n", if (spared2) "yes" else "NO"))
  cat(sprintf("   plant: a label with a mask drawn under it is spared: %s\n", if (spared3) "yes" else "NO"))
  quit(status = if (caught && spared1 && spared2 && spared3) 0 else 1)
}

.probe.files <- .probe.args[grepl("\\.R$", .probe.args)]
if (!length(.probe.files)) .probe.files <- Sys.glob("checks/fig_*.R")
for (.probe.f in .probe.files) {
  .ink$labels <- list(); .ink$marks <- list(); .ink$covers <- list(); .ink$name <- basename(.probe.f)
  .probe.dev <- tempfile(fileext = ".pdf"); pdf(.probe.dev, width = 8, height = 5)
  try(suppressWarnings(source(.probe.f, local = FALSE)), silent = TRUE)
  while (dev.cur() > 1) dev.off()
  .report(); unlink(.probe.dev)
}
if (length(.ink$hits)) {
  for (h in unique(.ink$hits)) cat("  ", h, "  <-- ISSUE\n")
  cat(sprintf("   %d label(s) printed over\n", length(unique(.ink$hits))))
  quit(status = 1)
}
cat(sprintf("   %d generators, no label printed over by the figure's own ink", length(.probe.files)))
if (length(.ink$ok)) cat(sprintf(" (%d accepted, with reasons in %s)", length(unique(.ink$ok)), .acc.f))
cat("\n")

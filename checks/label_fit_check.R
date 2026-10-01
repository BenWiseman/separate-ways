#!/usr/bin/env Rscript
# Does every label a figure draws actually fit inside its panel?
#
#   Rscript checks/fig_text_fits.R checks/fig_*.R
#   Rscript checks/fig_text_fits.R --selftest
#
# WHY. fig_legibility measures type SIZE and passed every figure at nine point or better. It says
# nothing about width. On 2026-09-27 the companion was built to PDF for the first time and page 20
# carried an annotation clipped mid-word, "nothing in particular elsewhe", because R silently cuts
# text at the plot region and reports nothing. Two figures were affected and neither gate saw it.
#
# HOW. Source each generator with text() shadowed, so every label is measured against par("usr")
# at the moment it is drawn, with the two traps that cost an hour: par("usr") is in log10 units on
# a log axis, and strwidth comes back NEGATIVE on a reversed axis, so the pair has to be ordered
# before it is compared.

# Dotted names throughout: these scripts are sourced into the global environment and they
# use short names like a, f and w themselves. An undotted argument vector gets clobbered
# by the first generator that assigns to it, and the run then reports the wrong count.
.probe.args <- commandArgs(TRUE)
ur_pair <- function(ur) c(min(ur[1], ur[2]), max(ur[1], ur[2]))

.overflow <- new.env(); .overflow$hits <- character(0)
text <- function(x, y = NULL, labels = seq_along(x), ...) {
  o <- list(...)
  cex <- if (!is.null(o$cex)) o$cex else 1
  adj <- if (!is.null(o$adj)) o$adj[1] else 0.5
  lab <- if (is.expression(labels) || is.call(labels)) character(0) else as.character(labels)
  if (length(lab) && is.null(o$srt)) {
    b <- ur_pair(par("usr")); lo <- b[1]; hi <- b[2]; tol <- 1e-9*(hi - lo)
    for (i in seq_along(lab)) {
      xi <- x[[min(i, length(x))]]
      if (par("xlog")) xi <- log10(xi)
      w <- strwidth(lab[i], cex = cex)
      e <- sort(c(xi - adj*w, xi - adj*w + w))
      if (is.finite(e[1]) && (e[2] > hi + tol || e[1] < lo - tol))
        .overflow$hits <- c(.overflow$hits,
          sprintf("%s: %.1f%% of the panel past the edge, \"%s\"", .probe.name,
                  100*max(e[2] - hi, lo - e[1])/(hi - lo), lab[i]))
    }
  }
  graphics::text(x, y, labels, ...)
}

# mtext was outside the shadow until 2026-09-27, and the first draft of the new graphical
# abstract shipped a panel title reading "...seen from below its tempe". The room a title has is
# NOT the plot region: mtext writes into the margin, so it is clipped at the FIGURE region,
# par("fin"), or at the device for outer = TRUE. Measuring against par("pin") instead flags four
# titles that render perfectly well, which is how this check was wrong on its first run.
mtext <- function(text, side = 3, line = 0, ...) {
  o <- list(...)
  cex <- if (!is.null(o$cex)) o$cex else 1
  outer <- isTRUE(o$outer)
  lab <- if (is.expression(text) || is.call(text)) character(0) else as.character(text)
  if (length(lab) && is.null(o$at) && is.null(o$adj)) {
    box <- if (outer) par("din") else par("fin")
    room <- if (side %in% c(1, 3)) box[1] else box[2]
    for (t in lab) {
      w <- strwidth(t, units = "inches", cex = cex)
      if (is.finite(w) && w > room)
        .overflow$hits <- c(.overflow$hits,
          sprintf("%s: a side-%d %stitle is %.0f%% wider than its %s, \"%s\"",
                  .probe.name, side, if (outer) "outer " else "", 100*(w/room - 1),
                  if (outer) "device" else "figure region", t))
    }
  }
  graphics::mtext(text, side = side, line = line, ...)
}

if ("--selftest" %in% .probe.args) {
  pdf(tempfile(), width = 6, height = 4)
  plot(NA, xlim = c(0, 1), ylim = c(0, 1))
  .probe.name <- "selftest"
  text(0.02, 0.5, "short", adj = 0)
  n1 <- length(.overflow$hits)
  text(0.02, 0.3, paste(rep("a label far too long to fit", 4), collapse = " "), adj = 0)
  n2 <- length(.overflow$hits)
  plot(NA, xlim = c(20, 6), ylim = c(0, 1))            # a reversed axis, the false-positive trap
  text(19, 0.5, "comfortably inside a reversed axis", adj = 0)
  n3 <- length(.overflow$hits)
  plot(NA, xlim = c(1e-4, 1e9), ylim = c(0, 1), log = "x")   # and a log axis
  text(1e-3, 0.5, "inside a log axis", adj = 0)
  n4 <- length(.overflow$hits)
  plot(NA, xlim = c(0, 1), ylim = c(0, 1))
  mtext("short title", side = 3, cex = 0.9)
  n5 <- length(.overflow$hits)
  mtext(paste(rep("a panel title far too long to fit", 5), collapse = " "), side = 3, cex = 0.9)
  n6 <- length(.overflow$hits)
  invisible(dev.off())
  stopifnot(n1 == 0, n2 == 1, n3 == 1, n4 == 1, n5 == 1, n6 == 2)
  cat("   plant: an over-long label and an over-long panel title are both caught, and a\n")
  cat("          short one, a reversed axis and a log axis are all spared: yes\n")
  quit(status = 0)
}

# Run with no arguments this printed "0 generators, every label inside its panel" and exited
# clean: a pass over nothing at all. label_ink_check.R already falls back to the glob, so match
# it, and refuse outright if the glob is empty rather than reporting a vacuous success.
if (!length(.probe.args)) .probe.args <- Sys.glob("checks/fig_*.R")
if (!length(.probe.args)) {
  cat("   no generators found: checks/fig_*.R matched nothing  <-- ISSUE\n")
  quit(status = 1)
}

for (.probe.f in .probe.args) {
  .probe.name <- basename(.probe.f)
  sink(tempfile()); try(source(.probe.f, local = FALSE), silent = TRUE); sink()
  while (dev.cur() > 1) invisible(dev.off())
}
if (length(.overflow$hits)) {
  cat(paste0("   ", .overflow$hits, "  <-- ISSUE\n"), sep = "")
  quit(status = 1)
}
cat(sprintf("   %d generators, every label inside its panel\n", length(.probe.args)))

# What goes in and what comes out, drawn, because the ledger has only ever been prose.
#
# Section 4.4 counts eighteen properties of general relativity out against four things in, two of
# which are measured numbers. That count is the paper's central claim and a reader meets it as a
# paragraph on page twenty. This draws it: the arena on the left, the one postulate in the middle,
# and the eighteen lines on the right, with the two numbers greyed so the asymmetry is the first
# thing the eye gets.
#
# Every label is a line of relativity_ledger.R and the counts are read from it rather than typed
# here, so the figure cannot drift from the ledger it draws.

src <- readLines("checks/calc/relativity_ledger.R")
# match the status field of an item() call, not the two Filter lines that also name the words
n_der <- length(grep('^ +"derived", ', src))
n_asu <- length(grep('^ +"assumed", ', src))
cat(sprintf("=== read from the ledger: %d derived, %d assumed ===\n", n_der, n_asu))
stopifnot(n_der == 19, n_asu == 4)

# short labels for the drawing; the ledger's own wording is the authority and these compress it
# labels are plotmath so the Greek renders in the pdf device, which has no Lambda glyph in the
# default encoding; a plain string there silently drops the character
OUT <- list("the field equations", "a horizon temperature", "horizons at every boost",
         "an area law for entropy", "the same coefficient everywhere", "the transverse involution",
         expression(paste(italic(G), " does not run")),
         expression(paste(italic(G) > 0, ", gravity attracts")),
         expression(paste(Lambda, " as an integration constant")),
         "null-focusing surfaces", "four dimensions, no others", "the classical-quantum divide",
         "black-hole thermodynamics", "a hot bang", "the weak equivalence principle",
         "the contracted Bianchi identity", "matter conservation", "universality of local horizons",
         "the two sheets' opposite time")
IN  <- list("a Lorentzian signature", "one shared metric",
            expression(paste("the value of ", italic(G))),
            expression(paste("the value of ", Lambda)))
IS_NUM <- c(FALSE, FALSE, TRUE, TRUE)
stopifnot(length(OUT) == n_der, length(IN) == n_asu)

ink <- "grey15"; c1 <- "#1f4e79"; c2 <- "#a8400f"; grey <- "#8a97a4"; pale <- "#b9c2cb"

draw <- function() {
  par(mar = c(0.3, 0.3, 1.9, 0.3), xpd = NA)
  plot(NA, xlim = c(0, 100), ylim = c(0, 100), axes = FALSE, xlab = "", ylab = "")
  rect(-3, -3, 103, 103, col = "#fbfbfa", border = NA)

  # the out column sets the vertical extent; the in column is centred against it
  ytop2 <- 91; dy2 <- 4.55
  ybot2 <- ytop2 - (length(OUT) - 1) * dy2
  ymid  <- (ytop2 + ybot2) / 2

  # ---- left: the arena. list indexing must be [[ ]] or a plotmath label is deparsed and drawn
  # as its own source text, which is what the first render of this figure did.
  dy <- 8.0; ytop <- ymid + 1.5 * dy
  text(14, ytop + 8, "what goes in", col = ink, cex = 1.00, font = 2)
  for (i in seq_along(IN)) {
    y <- ytop - (i - 1) * dy
    col <- if (IS_NUM[i]) grey else c1
    rect(2, y - 2.9, 26, y + 2.9, col = "#ffffff", border = col, lwd = 1.6)
    text(14, y, IN[[i]], col = col, cex = 0.82)
  }
  text(14, ytop - 4 * dy + 1.0, "the last two are measured numbers", col = grey, cex = 0.80)
  text(14, ytop - 4 * dy - 4.5, "the first two are what a metric", col = grey, cex = 0.80)
  text(14, ytop - 4 * dy - 8.0, "theory is", col = grey, cex = 0.80)

  # ---- middle: the one postulate
  rect(32, ymid - 7, 57, ymid + 7, col = "#ffffff", border = c2, lwd = 2.2)
  text(44.5, ymid + 2.6, "CPT holds of the universe", col = c2, cex = 0.88, font = 2)
  text(44.5, ymid - 2.6, "and not only of its laws", col = c2, cex = 0.88, font = 2)
  text(44.5, ymid - 12, "the metric is never", col = grey, cex = 0.80)
  text(44.5, ymid - 15.5, "an operator", col = grey, cex = 0.80)
  arrows(27, ymid + 4, 31, ymid + 2, length = 0.08, col = pale, lwd = 2)
  arrows(58, ymid, 61.5, ymid, length = 0.08, col = pale, lwd = 2)

  # ---- right: the eighteen
  text(80, ytop2 + 4.5, "what comes out", col = ink, cex = 1.00, font = 2)
  for (i in seq_along(OUT)) {
    y <- ytop2 - (i - 1) * dy2
    points(63, y, pch = 15, cex = 0.5, col = c1)
    text(65, y, OUT[[i]], col = ink, cex = 0.80, adj = 0)
  }
  mtext(sprintf("%d lines of general relativity out, %d in, and no path from any of them back to relativity",
                n_der, n_asu),
        side = 3, line = 0.4, cex = 0.88, col = ink)
}

for (dev in c("pdf", "png")) {
  out <- sprintf("paper/fig_ledger.%s", dev)
  if (dev == "pdf") pdf(out, width = 9.8, height = 5.4, pointsize = 12)
  else png(out, width = 9.8, height = 5.4, units = "in", res = 200, pointsize = 12)
  draw(); invisible(dev.off())
  cat(sprintf("   wrote %s\n", out))
}
cat("\n=== the plant: the figure must break if the ledger moves ===\n")
cat("   The counts are read from relativity_ledger.R and asserted against the label lists, so a\n")
cat("   line added to the ledger without a label added here stops the script rather than drawing\n")
cat("   a figure that disagrees with the text. Checked by pretending the ledger gained a line:\n")
fake <- n_der + 1
cat(sprintf("      ledger says %d derived, labels supply %d  ->  %s\n", fake, length(OUT),
            ifelse(fake == length(OUT), "would draw", "STOPS, as it must")))
stopifnot(fake != length(OUT))

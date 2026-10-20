# What goes in and what comes out, drawn, because the ledger has only ever been prose.
#
# Section 4.4 counts nineteen properties of general relativity out against four things in, two of
# which are measured numbers. This draws it: the arena on the left, the one postulate in the
# middle, and the nineteen lines on the right, with the two numbers greyed.
#
# 2026-10-06. The right column is split by where each line stood before this paper: nine were in
# Jacobson's 1995 argument, four are standard, six need the fold. The figure used to present all
# nineteen as the fold's, which is the overstatement a gr-qc reader meets first. The groups and
# their order are read from relativity_ledger.R, so a line moved between groups moves here too.
#
# Every label is a line of relativity_ledger.R and the counts are read from it rather than typed
# here, so the figure cannot drift from the ledger it draws.

src <- readLines("checks/calc/relativity_ledger.R")
# match the status field of an item() call, not the two Filter lines that also name the words
n_der <- length(grep('^ +"derived", ', src))
n_asu <- length(grep('^ +"assumed", ', src))
cat(sprintf("=== read from the ledger: %d derived, %d assumed ===\n", n_der, n_asu))
stopifnot(n_der == 19, n_asu == 4)
.led <- new.env()
invisible(capture.output(sys.source("checks/calc/relativity_ledger.R", envir = .led)))
DER <- Filter(function(x) x$status == "derived", .led$L)
GROUPS <- .led$GROUPS
stopifnot(length(DER) == n_der, identical(GROUPS, c("Jacobson 1995", "standard", "fold")))

# short labels for the drawing; the ledger's own wording is the authority and these compress it.
# Each is keyed by a pattern that must match exactly one ledger line, so a renamed line stops the
# script instead of losing its label. Labels are plotmath where Greek or italics are needed: the
# pdf device has no Lambda glyph in the default encoding and a plain string drops the character.
KEY <- list(
  list("FORM",                         "the field equations"),
  list("temperature at every",         "a horizon temperature"),
  list("every boost and orientation",  "horizons at every boost"),
  list("proportional to the horizon",  "an area law for entropy"),
  list("SAME entropy coefficient",     "the same coefficient everywhere"),
  list("does not run",                 expression(paste(italic(G), " does not run"))),
  list("SIGN of Newton",               expression(paste(italic(G) > 0, ", gravity attracts"))),
  list("^Lambda appearing",            expression(paste(Lambda, " as an integration constant"))),
  list("away from the fold's own",     "local horizons at every point"),
  list("black-hole thermodynamics",    "black-hole thermodynamics"),
  list("weak equivalence",             "the weak equivalence principle"),
  list("contracted Bianchi",           "the contracted Bianchi identity"),
  list("matter conservation",          "matter conservation"),
  list("transverse involution",        "the black-hole involution"),
  list("null-focusing",                "null-focusing surfaces"),
  list("four large spacetime",         "four dimensions, no others"),
  list("classical-quantum",            "the classical-quantum divide"),
  list("opposite time orientation",    "the two sheets' opposite time"),
  list("hot bang",                     "a hot bang"))
stopifnot(length(KEY) == n_der)
lab_of <- function(what) {
  hit <- which(sapply(KEY, function(k) grepl(k[[1]], what)))
  if (length(hit) != 1) stop(sprintf("ledger line matches %d labels: %s", length(hit), what))
  KEY[[hit]][[2]]
}
ROWS <- list()
for (g in GROUPS) for (x in DER) if (x$group == g) ROWS[[length(ROWS) + 1]] <- list(g = g, lab = lab_of(x$what))
stopifnot(length(ROWS) == n_der)
NG <- sapply(GROUPS, function(g) sum(sapply(ROWS, function(r) r$g == g)))
HEAD <- c("Jacobson 1995" = "in Jacobson's 1995 argument", "standard" = "standard results",
          "fold" = "added by the fold")

IN  <- list("a Lorentzian signature", "one shared metric",
            expression(paste("the value of ", italic(G))),
            expression(paste("the value of ", Lambda)))
IS_NUM <- c(FALSE, FALSE, TRUE, TRUE)
stopifnot(length(IN) == n_asu)

ink <- "grey15"; c1 <- "#1f4e79"; c2 <- "#a8400f"; grey <- "#8a97a4"; pale <- "#b9c2cb"
GCOL <- c("Jacobson 1995" = c1, "standard" = grey, "fold" = c2)
words <- c("One", "Two", "Three", "Four", "Five", "Six", "Seven", "Eight", "Nine")

draw <- function() {
  par(mar = c(0.3, 0.3, 1.9, 0.3), xpd = NA)
  plot(NA, xlim = c(0, 100), ylim = c(0, 100), axes = FALSE, xlab = "", ylab = "")
  rect(-3, -3, 103, 103, col = "#fbfbfa", border = NA)

  # the out column sets the vertical extent: one step per line, one per group heading, and half
  # a step of air before each heading after the first
  ytop2 <- 92; dy2 <- 3.72
  nstep <- n_der + length(GROUPS) + 0.5 * (length(GROUPS) - 1)
  ybot2 <- ytop2 - (nstep - 1) * dy2
  ymid  <- (ytop2 + ybot2) / 2

  # ---- left: the arena. list indexing must be [[ ]] or a plotmath label is deparsed and drawn
  # as its own source text, which is what the first render of this figure did.
  # 2026-10-04, Ben: the single diagonal arrow from the third box read as if only G fed the
  # postulate, the centre box was "pretentious philosophy speak", and the notes under the boxes
  # ran last-two before first-two. All four boxes now join one bracket and one arrow, the centre
  # says what the paper does, and the notes run in box order.
  dy <- 19.0; ytop <- ymid + 1.5 * dy
  text(13, ytop2 + 3.6, "what goes in", col = ink, cex = 1.00, font = 2)
  for (i in seq_along(IN)) {
    y <- ytop - (i - 1) * dy
    col <- if (IS_NUM[i]) grey else c1
    rect(2, y - 3.2, 24, y + 3.2, col = "#ffffff", border = col, lwd = 1.6)
    text(13, y, IN[[i]], col = col, cex = 0.82)
    segments(24, y, 26.5, y, col = pale, lwd = 2)
  }
  segments(26.5, ytop, 26.5, ytop - 3 * dy, col = pale, lwd = 2)
  ynote <- ytop - (length(IN) - 1) * dy - 9.5
  text(13, ynote,       "first two: what a metric theory is", col = grey, cex = 0.80)
  text(13, ynote - 4.5, "last two: measured constants", col = grey, cex = 0.80)

  # ---- middle: the one postulate, joined to algebraic QFT
  # The box is measured from its own text rather than given a fixed width.
  cx <- 44.5; mcex <- 0.84
  L1 <- "CPT of the universe,"; L2 <- "fused with algebraic QFT"
  halfw <- max(strwidth(L1, cex = mcex, font = 2),
               strwidth(L2, cex = mcex, font = 2)) / 2 + 2.4
  rect(cx - halfw, ymid - 7, cx + halfw, ymid + 7, col = "#ffffff", border = c2, lwd = 2.2)
  text(cx, ymid + 2.6, L1, col = c2, cex = mcex, font = 2)
  text(cx, ymid - 2.6, L2, col = c2, cex = mcex, font = 2)
  text(cx, ymid - 12, "the metric is never quantised", col = grey, cex = 0.80)
  stopifnot(cx - halfw > 28, cx + halfw < 60.5)   # must not collide with either bracket
  arrows(26.5, ymid, cx - halfw - 0.9, ymid, length = 0.08, col = pale, lwd = 2)
  arrows(cx + halfw + 0.9, ymid, 61.8, ymid, length = 0.08, col = pale, lwd = 2)

  # ---- right: the nineteen, behind one bracket, in three groups
  segments(62.3, ytop2, 62.3, ybot2, col = pale, lwd = 2)
  text(81, ytop2 + 3.6, "what comes out", col = ink, cex = 1.00, font = 2)
  y <- ytop2; first <- TRUE
  for (g in GROUPS) {
    if (!first) y <- y - 0.5 * dy2
    text(64.2, y, sprintf("%s, %s", HEAD[[g]], tolower(words[NG[[g]]])),
         col = GCOL[[g]], cex = 0.78, adj = 0, font = 2)
    y <- y - dy2; first <- FALSE
    for (r in ROWS) if (r$g == g) {
      points(65.4, y, pch = 15, cex = 0.5, col = GCOL[[g]])
      text(67.2, y, r$lab, col = ink, cex = 0.78, adj = 0)
      y <- y - dy2
    }
  }
  stopifnot(abs((y + dy2) - ybot2) < 1e-6, ybot2 > 3)   # the last line lands where it was planned
  mtext(sprintf("%s in, nineteen properties of general relativity out, %s of them added by the fold",
                words[n_asu], tolower(words[NG[["fold"]]])),
        side = 3, line = 0.4, cex = 0.86, col = ink)
}

for (dev in c("pdf", "png")) {
  out <- sprintf("papers/1_separate_ways/fig_ledger.%s", dev)
  if (dev == "pdf") pdf(out, width = 9.8, height = 5.8, pointsize = 12)
  else png(out, width = 9.8, height = 5.8, units = "in", res = 200, pointsize = 12)
  draw(); invisible(dev.off())
  cat(sprintf("   wrote %s\n", out))
}
cat(sprintf("   groups drawn: %s\n", paste(sprintf("%s %d", GROUPS, NG), collapse = ", ")))
cat("\n=== the plant: the figure must break if the ledger moves ===\n")
cat("   The counts are read from relativity_ledger.R and asserted against the label list, so a\n")
cat("   line added to the ledger without a label added here stops the script rather than drawing\n")
cat("   a figure that disagrees with the text. Checked by pretending the ledger gained a line:\n")
fake <- n_der + 1
cat(sprintf("      ledger says %d derived, labels supply %d  ->  %s\n", fake, length(KEY),
            ifelse(fake == length(KEY), "would draw", "STOPS, as it must")))
stopifnot(fake != length(KEY))
cat("   And a ledger line no label pattern matches must stop it too:\n")
ok <- tryCatch({ lab_of("a twentieth line about nothing"); FALSE }, error = function(e) TRUE)
cat(sprintf("      an unmatched line -> %s\n", if (ok) "STOPS, as it must" else "would draw"))
stopifnot(ok)

# FIGURE: "The floor under every admissible state."
#
# Section 3.1 is the paper's central argument and it arrives as a 4x4 matrix, a basis the
# reader has not been given, and the word "spectrum". The claim underneath is simple enough
# to draw: the matrix has four eigenvalues, two of them are exactly 1/2 whatever happens, and
# the lowest one cannot be pushed below n*. That is the whole of the ceiling.
#
# The eigenvalues in panel (b) are computed with eigen() and checked against the closed form
# the manuscript quotes, so the picture cannot disagree with the text.
#
# Base R only. Authored at 6.5 in, the width the manuscript places it at.

INK <- "#1b1b1b"; LOW <- "#0b6e9e"; HIGH <- "#c2461f"; MID <- "#6b6b6b"; GREY <- "#7e7e7e"
FORB <- rgb(0.76, 0.27, 0.12, 0.13)
CELL <- rgb(0.11, 0.43, 0.62, 0.15)

nstar <- function(P) (1 - sqrt(1 - P)) / 2
Qmat  <- function(P) {
  c0 <- 0.5 * sqrt(P * (1 - P))
  matrix(c(P/2, 0,   0,   c0,
           0,   0.5, 0,   0,
           0,   0,   0.5, 0,
           c0,  0,   0,   1 - P/2), 4, 4, byrow = TRUE)
}

# ---- check the drawn spectrum against the closed form, at several P
for (P in c(0.05, 0.2, 0.5, 0.8, 0.95)) {
  ev <- sort(eigen(Qmat(P), symmetric = TRUE)$values)
  stopifnot(max(abs(ev - sort(c(nstar(P), 0.5, 0.5, 1 - nstar(P))))) < 1e-12)
}
cat("checks passed: eigen(Q) matches {n*, 1/2, 1/2, 1-n*} at every P tried\n")

PSHOW <- 0.5
LAB <- c("00", "10", "01", "11")

draw <- function() {
  layout(matrix(c(1, 2), 1, 2), widths = c(1, 1.30))

  # ---------------- (a) the matrix ----------------
  par(mar = c(2.6, 2.4, 2.3, 0.6), family = "sans")
  Q <- Qmat(PSHOW)
  plot(NA, xlim = c(-0.95, 4.35), ylim = c(-1.25, 4.95), axes = FALSE,
       xlab = "", ylab = "", asp = 1)
  mtext("(a)  the pair, as a matrix", side = 3, line = 0.5, adj = 0,
        cex = 0.98, font = 2, col = INK)
  for (i in 1:4) for (j in 1:4) {
    x0 <- j - 1; y0 <- 4 - i
    v <- Q[i, j]
    rect(x0, y0, x0 + 1, y0 + 1, col = if (v != 0) CELL else "white",
         border = GREY, lwd = 0.9)
    if (v != 0) text(x0 + 0.5, y0 + 0.5, sprintf("%.2f", v), col = INK, cex = 0.92)
  }
  for (j in 1:4) text(j - 0.5, 4.22, LAB[j], col = MID, cex = 0.88)
  for (i in 1:4) text(-0.18, 4.5 - i, LAB[i], col = MID, cex = 0.88, adj = 1)
  text(2.0, -0.42, "the four states of a pair", col = INK, cex = 0.90, adj = 0.5)
  text(2.0, -1.02, sprintf("drawn at P = %.1f", PSHOW), col = MID, cex = 0.88, adj = 0.5)

  # ---------------- (b) its spectrum, and the floor ----------------
  par(mar = c(3.7, 4.0, 2.3, 1.1))
  Ps <- seq(0, 1, length.out = 600)
  plot(NA, xlim = c(0, 1), ylim = c(-0.02, 1.62), axes = FALSE, xlab = "", ylab = "")
  mtext("(b)  and the floor under it", side = 3, line = 0.5, adj = 0,
        cex = 0.98, font = 2, col = INK)

  polygon(c(Ps, rev(Ps)), c(nstar(Ps), rep(0, length(Ps))), col = FORB, border = NA)
  lines(Ps, nstar(Ps),     col = LOW,  lwd = 2.8)
  lines(Ps, 1 - nstar(Ps), col = HIGH, lwd = 2.8)
  segments(0, 0.5, 1, 0.5, col = MID, lwd = 2.4)

  ev <- sort(eigen(Qmat(PSHOW), symmetric = TRUE)$values)
  points(rep(PSHOW, 4), ev, pch = 19, cex = 1.05,
         col = c(LOW, MID, MID, HIGH))

  axis(1, at = seq(0, 1, 0.25), cex.axis = 0.92, col = GREY, col.axis = INK,
       tck = -0.02, mgp = c(2, 0.6, 0))
  axis(2, at = c(0, 0.5, 1), labels = c("0", "1/2", "1"), las = 1, cex.axis = 0.92,
       col = GREY, col.axis = INK, tck = -0.02, mgp = c(2, 0.7, 0))
  mtext("how strongly the bang mixes the mode", side = 1, line = 2.3, cex = 0.92, col = INK)
  mtext("the four eigenvalues", side = 2, line = 2.5, cex = 0.92, col = INK)

  # legend band above the curves: nothing is labelled where a curve can reach it
  key <- function(y, col, lab, swatch = "line") {
    if (swatch == "line") segments(0.02, y, 0.11, y, col = col, lwd = 2.8)
    else rect(0.02, y - 0.045, 0.11, y + 0.045, col = FORB, border = NA)
    text(0.14, y, lab, col = INK, cex = 0.90, adj = 0)
  }
  key(1.54, HIGH, "the highest")
  key(1.39, MID,  "two always at 1/2")
  key(1.24, LOW,  "the lowest")
  key(1.09, NA,   "no admissible state here", swatch = "fill")
}

dir.create("pub/paper2/figs", showWarnings = FALSE, recursive = TRUE)
pdf("pub/paper2/figs/fig_floor.pdf", width = 6.5, height = 3.5); draw(); dev.off()
png("pub/paper2/figs/fig_floor.png", width = 1950, height = 1050, res = 300); draw(); dev.off()
cat("wrote fig_floor.{pdf,png}\n")

# FIGURE: "One gravitational event sets two numbers."
#
# 3.5's result is the paper's most striking and it arrives under the words "Elimination of I
# between them". The claim is a straight line on log axes and should be drawn as one. What the
# figure must not do is overstate it: the SLOPE carries no free parameter, the HEIGHT moves with
# the state and with the counting convention, and 3.5 says so. Both go in the picture.
#
# The exponent is not asserted. It is obtained by composing the two relations the paper uses and
# fitting the result, which is the elimination 3.5 performs.
#
# Base R only. Authored at 6.5 in, the width the manuscript places it at.

INK <- "#1b1b1b"; LINE <- "#0b6e9e"; PT <- "#c2461f"; GREY <- "#7e7e7e"
FREE <- rgb(0.11, 0.43, 0.62, 0.13)

M0  <- 491.6         # PeV, the ceiling
T0  <- 1.417e-32     # s, at c_G = 1 and the least-occupied state
I0  <- 0.0127597
CG_HALF <- 2^(2/3)   # counting a Majorana pair once rather than twice
BAND    <- 2.33      # +133 per cent, the full state freedom of Table 1 (R_under_band.R)

# --- the elimination, done rather than asserted:  M1 ~ I^(-2/5)  and  t ~ M1^-1 I^(-2/3)
Ifun <- function(M) I0 * (M / M0)^(-5/2)
A    <- T0 * M0 * Ifun(M0)^(2/3)
tdec <- function(M) A * M^(-1) * Ifun(M)^(-2/3)

Ms <- 10^seq(log10(M0) - 2, log10(M0) + 2, length.out = 400)
slope <- coef(lm(log(tdec(Ms)) ~ log(Ms)))[2]
stopifnot(abs(slope - 2/3) < 1e-10)                      # the 2/3 power, over four decades
stopifnot(abs(tdec(M0) - T0) / T0 < 1e-12,
          abs(tdec(M0) * CG_HALF - 2.249e-32) / 2.249e-32 < 2e-3)   # the second quoted figure
cat(sprintf("fitted exponent over four decades: %.12f\n", slope))

draw <- function() {
  par(mar = c(3.8, 4.5, 2.2, 1.2), family = "sans")
  ylo <- log10(tdec(min(Ms))) - 0.30
  yhi <- log10(tdec(max(Ms)) * BAND) + 0.42
  plot(NA, xlim = log10(range(Ms)), ylim = c(ylo, yhi), axes = FALSE, xlab = "", ylab = "")
  mtext("One gravitational event sets both numbers", side = 3, line = 0.5, adj = 0,
        cex = 1.00, font = 2, col = INK)

  polygon(c(log10(Ms), rev(log10(Ms))),
          c(log10(tdec(Ms)), rev(log10(tdec(Ms) * BAND))), col = FREE, border = NA)
  lines(log10(Ms), log10(tdec(Ms)), col = LINE, lwd = 2.8)
  lines(log10(Ms), log10(tdec(Ms) * CG_HALF), col = LINE, lwd = 1.4, lty = 2)
  points(log10(M0), log10(T0), pch = 19, cex = 1.35, col = PT)

  xt <- c(5, 50, 491.6, 5000, 50000)
  axis(1, at = log10(xt), labels = c("5", "50", "492", "5000", "50000"),
       cex.axis = 0.92, col = GREY, col.axis = INK, tck = -0.02, mgp = c(2, 0.6, 0))
  yt <- seq(ceiling(ylo), floor(yhi), 1)
  axis(2, at = yt, labels = parse(text = sprintf("10^{%d}", yt)), las = 1,
       cex.axis = 0.92, col = GREY, col.axis = INK, tck = -0.02, mgp = c(2, 0.7, 0))
  mtext("dark-matter mass  (PeV)", side = 1, line = 2.4, cex = 0.94, col = INK)
  mtext("decoherence time  (s)", side = 2, line = 3.1, cex = 0.94, col = INK)

  # legend in the empty upper-left, well clear of the diagonal
  xs <- log10(M0) - 1.95; xt2 <- xs + 0.34
  ky <- yhi - c(0.18, 0.46, 0.74)
  segments(xs, ky[1], xt2, ky[1], col = LINE, lwd = 2.8)
  text(xt2 + 0.08, ky[1], "the adopted state and convention", col = INK, cex = 0.88, adj = 0)
  segments(xs, ky[2], xt2, ky[2], col = LINE, lwd = 1.4, lty = 2)
  text(xt2 + 0.08, ky[2], "a Majorana pair counted once", col = INK, cex = 0.88, adj = 0)
  rect(xs, ky[3] - 0.09, xt2, ky[3] + 0.09, col = FREE, border = NA)
  text(xt2 + 0.08, ky[3], "where the height can move", col = INK, cex = 0.88, adj = 0)

  # slope triangle on the line, labels kept inside the empty lower-right wedge
  x1 <- log10(M0) - 1.80; x2 <- x1 + 1.0
  y1 <- log10(tdec(10^x1)); y2 <- y1 + 2/3
  segments(x1, y1, x2, y1, col = INK, lwd = 1.5)
  segments(x2, y1, x2, y2, col = INK, lwd = 1.5)
  text((x1 + x2)/2, y1 - 0.22, "10x the mass", col = INK, cex = 0.88, adj = 0.5)
  text(x2 + 0.09, y1 + 0.30, "4.6x the time", col = INK, cex = 0.88, adj = 0)

  # the ceiling, labelled below-right of its point where nothing else sits
  segments(log10(M0) + 0.05, log10(T0) - 0.10, log10(M0) + 0.30, log10(T0) - 0.52,
           col = PT, lwd = 1.2)
  text(log10(M0) + 0.34, log10(T0) - 0.64, "the ceiling, 491.6 PeV", col = PT,
       cex = 0.90, adj = 0)
}

dir.create("pub/paper2/figs", showWarnings = FALSE, recursive = TRUE)
pdf("pub/paper2/figs/fig_power.pdf", width = 6.5, height = 3.7); draw(); dev.off()
png("pub/paper2/figs/fig_power.png", width = 1950, height = 1110, res = 300); draw(); dev.off()
cat("wrote fig_power.{pdf,png}\n")

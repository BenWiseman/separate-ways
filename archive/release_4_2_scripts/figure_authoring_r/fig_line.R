# FIGURE: "The line, the event, and the test that would settle it."
#
# Section 3.2 is the most testable thing in the paper and had no picture. It carries two
# statements a reader has to hold at once: the endpoint is sharp and the one measured event is
# not, so the event decides nothing; and a hundred events would decide it. Draw both.
#
# Panel (b) PARSES Table 2 out of the manuscript rather than repeating its numbers, so the
# figure cannot drift from the table the way hardcoded values have before.
#
# Base R only. Authored at 6.5 in, the width the manuscript places it at.

INK <- "#1b1b1b"; END <- "#c2461f"; EVT <- "#0b6e9e"; GREY <- "#7e7e7e"
ABOVE <- rgb(0.76, 0.27, 0.12, 0.17); BELOW <- rgb(0.11, 0.43, 0.62, 0.13)

E0 <- 245.8                       # the two-body endpoint, PeV
m <- 220; lo <- 110; hi <- 790    # KM3NeT reconstructed median and 68% interval [35]
s1 <- log(m/lo); s2 <- log(hi/m)

# Split-lognormal with each side carrying half the mass, which is the form
# separate_ways/tangents/stats/km3net_posterior.R uses and the form 3.2's number comes from.
# It is deliberately NOT the continuous-density split-normal: weighting the two halves to make
# the density continuous at the median moves 13.8 percentage points of mass across the endpoint
# and would have made the shaded area disagree with its own label. Both forms reproduce the
# published 68 per cent interval, so only the tail check below separates them.
dens <- function(E) {
  x <- log(E/m)
  ifelse(x < 0, dnorm(x, 0, s1), dnorm(x, 0, s2)) / E
}
p_above <- 1 - pnorm(log(E0/m)/s2)
stopifnot(abs(p_above - 0.465) < 0.005)          # the 46 per cent 3.2 quotes

# The curve and the annotation must be the same distribution. This is the check that was
# missing: it compares the AREA actually drawn against the number actually printed.
# The equal-mass form is discontinuous at the median by construction, so every integral is
# split there; quadrature run straight across the jump is off by 4e-6 rather than 7e-9.
.mass <- function(a, b) {
  if (a < m && b > m) integrate(dens, a, m)$value + integrate(dens, m, b)$value
  else                integrate(dens, a, b)$value
}
stopifnot(abs(.mass(1e-6, Inf) - 1)       < 1e-7,   # normalised
          abs(.mass(E0,   Inf) - p_above) < 1e-7,   # the area drawn IS the number printed
          abs(.mass(lo,    hi) - 0.6827)  < 1e-3)   # reproduces the published 68% interval

# ---- Table 2, read from the manuscript so the two cannot disagree
md <- readLines("pub/paper2/PAPER2_v3.md", warn = FALSE)
rows <- grep("^\\| *[0-9]+ *\\| *[0-9]+ *\\| *0?1?\\.[0-9]{2} *\\|", md, value = TRUE)
stopifnot(length(rows) == 4)
tab <- do.call(rbind, lapply(rows, function(r) {
  as.numeric(trimws(strsplit(r, "|", fixed = TRUE)[[1]][-1]))
}))
colnames(tab) <- c("N", "thresh", "r2", "r3", "r5")
cat("Table 2 as read from the manuscript:\n"); print(tab)
stopifnot(nrow(tab) == 4, tab[2, "N"] == 100, abs(tab[2, "r2"] - 0.88) < 1e-9)

draw <- function() {
  layout(matrix(c(1, 2), 1, 2), widths = c(1.15, 1.05))

  # ---------------- (a) the endpoint against the one event we have ----------------
  par(mar = c(3.7, 1.0, 2.3, 1.0), family = "sans")
  Es <- exp(seq(log(55), log(1700), length.out = 900))
  d <- dens(Es); dmax <- max(d)
  plot(NA, xlim = log(c(55, 1700)), ylim = c(-0.10, 1.58), axes = FALSE, xlab = "", ylab = "")
  mtext("(a)  the one event decides nothing", side = 3, line = 0.5, adj = 0,
        cex = 0.98, font = 2, col = INK)

  up <- Es >= E0; dn <- Es <= E0
  polygon(log(c(Es[dn], rev(Es[dn]))), c(d[dn], rep(0, sum(dn)))/dmax, col = BELOW, border = NA)
  polygon(log(c(Es[up], rev(Es[up]))), c(d[up], rep(0, sum(up)))/dmax, col = ABOVE, border = NA)
  lines(log(Es), d/dmax, col = EVT, lwd = 2.3)
  segments(log(E0), 0, log(E0), 1.18, col = END, lwd = 2.6)
  # 2026-09-23: the equal-mass split form steps at the median by construction, and with the
  # endpoint line only 26 PeV to its right the step read as a registration error rather than as a
  # feature. A light guide ties the step to the median marker so the reader sees it is deliberate.
  segments(log(m), 0, log(m), dens(m) / dmax, col = adjustcolor(EVT, 0.45), lwd = 1.1, lty = 3)
  points(log(m), 0, pch = 17, cex = 1.15, col = EVT)

  at <- c(100, 250, 500, 1000)
  axis(1, at = log(at), labels = at, cex.axis = 0.92, col = GREY, col.axis = INK,
       tck = -0.02, mgp = c(2, 0.6, 0))
  mtext("neutrino energy  (PeV)", side = 1, line = 2.3, cex = 0.92, col = INK)

  text(log(E0) + 0.06, 1.50, "the endpoint,", col = END, cex = 0.90, adj = 0)
  text(log(E0) + 0.06, 1.36, "245.8 PeV", col = END, cex = 0.90, adj = 0)
  text(log(E0) - 0.08, 1.50, "the event,", col = EVT, cex = 0.90, adj = 1)
  text(log(E0) - 0.08, 1.36, "as reconstructed", col = EVT, cex = 0.90, adj = 1)
  text(log(700), 0.46, sprintf("%.1f%% of it", 100 * p_above), col = END, cex = 0.90, adj = 0.5)
  text(log(700), 0.32, "sits above", col = END, cex = 0.90, adj = 0.5)
  text(log(m), -0.09, "median", col = EVT, cex = 0.88, adj = 0.5)

  # ---------------- (b) what a population would settle ----------------
  par(mar = c(3.7, 3.5, 2.3, 1.1))
  plot(NA, xlim = log(c(26, 1300)), ylim = c(0.52, 1.30), axes = FALSE, xlab = "", ylab = "")
  mtext("(b)  a hundred would decide", side = 3, line = 0.5, adj = 0,
        cex = 0.98, font = 2, col = INK)
  COLS <- c(END, "#8a6d1f", EVT); KEYS <- c("r2", "r3", "r5")
  LABS <- c("twice as high", "three times", "five times")
  for (k in 3:1) lines(log(tab[, "N"]), tab[, KEYS[k]], col = COLS[k], lwd = 2.5,
                       type = "b", pch = 19, cex = 0.85)

  axis(1, at = log(tab[, "N"]), labels = tab[, "N"], cex.axis = 0.92, col = GREY,
       col.axis = INK, tck = -0.02, mgp = c(2, 0.6, 0))
  axis(2, at = seq(0.6, 1.0, 0.2), las = 1, cex.axis = 0.92, col = GREY,
       col.axis = INK, tck = -0.02, mgp = c(2, 0.7, 0))
  mtext("events assigned to decay", side = 1, line = 2.3, cex = 0.92, col = INK)
  mtext("chance the test fires", side = 2, line = 2.5, cex = 0.92, col = INK)

  for (k in 1:3) {
    y <- 1.27 - (k - 1) * 0.075
    segments(log(30), y, log(52), y, col = COLS[k], lwd = 2.5)
    text(log(58), y, LABS[k], col = INK, cex = 0.88, adj = 0)
  }
}

dir.create("pub/paper2/figs", showWarnings = FALSE, recursive = TRUE)
pdf("pub/paper2/figs/fig_line.pdf", width = 6.5, height = 3.5); draw(); dev.off()
png("pub/paper2/figs/fig_line.png", width = 1950, height = 1050, res = 300); draw(); dev.off()
cat(sprintf("wrote fig_line.{pdf,png}   P(E > endpoint) = %.3f\n", p_above))

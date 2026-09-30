#!/usr/bin/env Rscript
# fig_overlap_threshold.R -- Section 2: when the two time directions stop being one state.
#
# The cosmology paper's 4.1 gives the overlap of the two branch states on closed de Sitter as
#     log |<E-|E+>| = -(N_f/4) sum_{n>=2} d_n log[1 + A^4(A^2-1) / (n^2 (n^2-1)^2)],   A = aH,
# and the whole question is which d_n. Counting every oscillator, d_n = n^2. Counting only the
# P_perp-invariant harmonics, d_n is the sum of (2l+1) over EVEN l, because P_perp acts as
# (-1)^l and so SPLITS each level rather than signing it.
#
# Section 2.1 settles it: the branch states are Gaussian in each mode coordinate, and a Gaussian
# in q^2 cannot tell +q from -q, so the projector acts as the identity on every mode whatever
# the parity of its harmonic. Restricting the physical states removes no oscillator, and the
# ungauged count is the one that applies. An earlier draft used the gauged one and put the
# threshold 19 per cent later, and that delay is withdrawn.
# Arithmetic: separate_ways/tangents/info/wdw_weights.R

dn_ungauged <- function(n) n^2
dn_gauged   <- function(n) { l <- 0:(n-1); sum((2*l+1)[l %% 2 == 0]) }
# d_n does not depend on A, so build both vectors once. Recomputing them inside the sum made
# this script take minutes instead of a second.
NMAX <- 4000
NN   <- 2:NMAX
DEN  <- NN^2 * (NN^2 - 1)^2
D_UN <- sapply(NN, dn_ungauged)
D_GA <- sapply(NN, dn_gauged)
logov <- function(A, d, Nf = 1) -(Nf/4) * sum(d * log(1 + A^4*(A^2-1)/DEN))
Aun <- uniroot(function(A) logov(A, D_UN) + 1, c(1.2, 3.5), tol = 1e-12)$root
Aga <- uniroot(function(A) logov(A, D_GA)   + 1, c(1.2, 4.5), tol = 1e-12)$root
cat("=== the decoherence threshold, both counts ===\n\n")
cat(sprintf("  ungauged (every oscillator):   aH = %.5f,  Ht = acosh(aH) = %.5f\n", Aun, acosh(Aun)))
cat(sprintf("  gauged   (invariant harmonics): aH = %.5f,  Ht = acosh(aH) = %.5f\n", Aga, acosh(Aga)))
cat(sprintf("  the gauged reading is %.1f per cent later in aH and %.1f per cent later in Ht,\n",
            100*(Aga/Aun - 1), 100*(acosh(Aga)/acosh(Aun) - 1)))
cat("  and it is the one withdrawn. The paper quotes the Ht figure, which is the proper time.\n")
stopifnot(abs(Aun - 1.95376) < 5e-6, abs(Aga - 2.43912) < 5e-6,
          abs(acosh(Aun) - 1.28984) < 1e-4, abs(acosh(Aga) - 1.53984) < 1e-4)
cat("\n  at n = 2 the four harmonics split one invariant and three not:")
cat(sprintf(" %d against %d\n", dn_gauged(2), dn_ungauged(2)))
stopifnot(dn_gauged(2) == 1, dn_ungauged(2) == 4)

ink <- "grey15"; c1 <- "#1f4e79"; c2 <- "#8a97a4"; mark <- "#a8400f"
draw <- function() {
  par(mar = c(3.9, 4.8, 2.2, 1.2), mgp = c(2.6, 0.7, 0), xpd = FALSE)
  A <- seq(1.0, 3.0, by = 0.004)
  un <- sapply(A, function(a) exp(logov(a, D_UN)))
  ga <- sapply(A, function(a) exp(logov(a, D_GA)))
  plot(NA, xlim = c(1, 3), ylim = c(0, 1.04), axes = FALSE,
       xlab = expression(paste("scale factor in Hubble units  ", italic(aH))),
       ylab = "overlap between the two branch states")
  rect(1, 0, 3, 1.04, col = "#fbfbfa", border = NA)
  axis(1, at = seq(1, 3, 0.5)); axis(2, at = seq(0, 1, 0.25), las = 1)
  segments(1, exp(-1), 3, exp(-1), col = mark, lwd = 1.6)
  text(1.03, exp(-1) + 0.048, expression(paste("one e-fold of overlap lost, ", italic(e)^-1)),
       col = mark, cex = 0.82, adj = 0)
  lines(A, ga, col = c2, lwd = 2.4, lty = 2)
  lines(A, un, col = c1, lwd = 2.8)
  points(Aun, exp(-1), pch = 19, cex = 1.3, col = c1)
  points(Aga, exp(-1), pch = 19, cex = 1.2, col = c2)
  segments(Aun, 0, Aun, exp(-1), col = c1, lty = 3); segments(Aga, 0, Aga, exp(-1), col = c2, lty = 3)
  text(Aun - 0.04, 0.115, "1.95376", col = c1, cex = 0.82, adj = 1)
  text(Aun - 0.04, 0.060, expression(paste(italic(Ht), " = 1.28984")), col = c1, cex = 0.82, adj = 1)
  text(Aga + 0.04, 0.115, "2.43912", col = c2, cex = 0.82, adj = 0)
  text(Aga + 0.04, 0.060, expression(paste(italic(Ht), " = 1.53984")), col = c2, cex = 0.82, adj = 0)
  # placed in the upper right, which both curves have left by aH = 2.2
  segments(1.80, 0.955, 1.96, 0.955, col = c1, lwd = 2.8)
  text(2.00, 0.955, "every oscillator counted, as 2.1 argues", col = c1, cex = 0.82, adj = 0)
  segments(1.80, 0.875, 1.96, 0.875, col = c2, lwd = 2.4, lty = 2)
  text(2.00, 0.875, "invariant harmonics only, withdrawn", col = c2, cex = 0.82, adj = 0)
  mtext("when the two time directions stop being one state", side = 3, line = 0.7,
        cex = 0.88, col = ink)
}
for (f in c("papers/2_over_the_horizon/fig_companion_overlap.pdf", "papers/2_over_the_horizon/fig_companion_overlap.png")) {
  if (grepl("pdf$", f)) pdf(f, width = 7.4, height = 4.4) else png(f, width = 1150, height = 685, res = 150)
  draw(); invisible(dev.off()); cat("  wrote", f, "\n")
}

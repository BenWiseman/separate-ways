#!/usr/bin/env Rscript
# fig_spin_branch.R -- Section 6: the fold-invariant class is a curve, and it has an exact end.
#
# Kerr-de Sitter has Delta_r = (r^2+a^2)(1 - Lam r^2/3) - 2 M r. Everything dimensionless:
# write r = M rho, a = M alpha, y = Lam M^2, so
#     Delta_r / M^2 = -(y/3) rho^4 + (1 - y alpha^2 / 3) rho^2 - 2 rho + alpha^2 .
# A fold-invariant hole needs the black-hole and cosmological horizons to coincide, which is
# a DOUBLE root of Delta_r. For each alpha that fixes y, so the class is a curve in the
# (a/M, 9 Lam M^2) plane rather than the single point A.15 finds at zero spin.
#
# The curve ends where the double root becomes triple. Matching
#     -(Lam/3)(r - r0)^3 (r + 3 r0)
# term by term gives Lam r0^2 = 2 sqrt(3) - 3 =: u, then M = (4/3) u r0 and a = sqrt(u) r0, so
#     a/M   = 3 / (4 sqrt(u)) = 1.1009173688      (the same in every de Sitter background)
#     9LamM^2 = 16 u^3
# Both are checked below against the numerical branch.
# Arithmetic and the isometry check: checks/calc/kerr_de_sitter_nariai.R

u     <- 2*sqrt(3) - 3
a_end <- 3 / (4 * sqrt(u))
y_end <- 16 * u^3                      # this is 9 Lam M^2, not Lam M^2
stopifnot(abs(a_end - 1.1009173688) < 1e-9)

# Solving the double-root condition by hunting roots of Delta_r fails near the end of the
# branch, where three of them coincide and polyroot loses two thirds of its digits. Eliminate
# instead. With A = 1 - y alpha^2/3, Delta = 0 and Delta' = 0 give A = (1 + (2y/3) rho^3)/rho
# from the second; substituting into the first collapses it to y = 3(rho - alpha^2)/rho^4, and
# putting both back into the definition of A leaves ONE quartic in rho alone:
#
#       rho^4 - 3 rho^3 + 2 alpha^2 rho^2 - alpha^2 rho + alpha^4 = 0 .
#
# At alpha = 0 that is rho^3(rho - 3), so rho = 3 and y = 1/9: Nariai, 9 Lam M^2 = 1. The branch
# is this quartic's largest real root, and it ENDS where that root disappears.

branch <- function(al) {
  cf <- c(al^4, -al^2, 2*al^2, -3, 1)          # constant first, for polyroot
  z  <- polyroot(cf); rr <- Re(z[abs(Im(z)) < 1e-6 * pmax(1, abs(Re(z)))])
  rr <- rr[rr > 0]
  if (!length(rr)) return(c(NA, NA))
  rho <- max(rr)
  y   <- 3 * (rho - al^2) / rho^4
  if (!is.finite(y) || y <= 0) return(c(NA, NA))
  c(9 * y, rho)
}

als <- c(seq(0, 1.09, length.out = 150), seq(1.0905, a_end, length.out = 90))
br  <- t(sapply(als, branch))
ys  <- br[,1]; rhos <- br[,2]
good <- is.finite(ys)
cat("=== the fold-invariant branch in Kerr-de Sitter ===\n\n")
cat(sprintf("  at a/M = 0        9 Lam M^2 = %.6f   (Nariai, exactly 1)\n", ys[1]))
last <- max(which(good))
cat(sprintf("  at a/M = %.6f  9 Lam M^2 = %.6f\n", als[last], ys[last]))
cat(sprintf("  closed form endpoint: a/M = %.10f, 9 Lam M^2 = 16 u^3 = %.6f\n", a_end, y_end))
cat(sprintf("  numerical endpoint agrees to %.2e in 9 Lam M^2\n", abs(ys[last] - y_end)))
stopifnot(abs(ys[1] - 1) < 1e-8, abs(ys[last] - y_end) < 1e-6, abs(als[last] - a_end) < 1e-9)

# The check must be able to fail: past the endpoint no double root exists at all.
cat(sprintf("\n  past the end, at a/M = %.4f: branch returns %s\n", a_end + 0.02,
            ifelse(is.na(branch(a_end + 0.02)[1]), "no merged pair, as it must", "SOMETHING")))
stopifnot(is.na(branch(a_end + 0.05)[1]))
cat(sprintf("  and the Kerr bound a/M = 1 is crossed at 9 Lam M^2 = %.4f, which a positive\n",
            approx(als[good], ys[good], 1)$y))
cat("  Lambda permits: extremal Kerr is not the limit here.\n")

ink <- "grey15"; c1 <- "#1f4e79"; mark <- "#a8400f"; grey <- "#8a97a4"
draw <- function() {
  par(mar = c(3.9, 4.6, 2.2, 1.4), mgp = c(2.5, 0.7, 0), xpd = FALSE)
  plot(NA, xlim = c(0, 1.22), ylim = c(0.90, 1.70), axes = FALSE,
       xlab = expression(paste("spin  ", italic(a), "/", italic(M))),
       ylab = expression(paste("9 ", Lambda, italic(M)^2)))
  rect(0, 0.90, 1.22, 1.70, col = "#fbfbfa", border = NA)
  axis(1, at = seq(0, 1.2, 0.2)); axis(2, at = seq(1.0, 1.6, 0.2), las = 1)

  rect(1, 0.90, 1.22, 1.70, col = "#f3f0ec", border = NA)
  segments(1, 0.90, 1, 1.70, col = grey, lty = 2)
  text(1.015, 1.075, "past extremal", col = grey, cex = 0.80, adj = 0)
  text(1.015, 1.030, "Kerr, which a", col = grey, cex = 0.80, adj = 0)
  text(1.015, 0.985, expression(paste("positive ", Lambda)), col = grey, cex = 0.80, adj = 0)
  text(1.015, 0.940, "permits", col = grey, cex = 0.80, adj = 0)

  lines(als[good], ys[good], col = c1, lwd = 2.8)

  points(a_end, y_end, pch = 19, cex = 1.45, col = mark)
  text(0.06, 1.665, expression(paste("the branch ends at  ", italic(a), "/", italic(M), " = 3/(4",
       sqrt(2*sqrt(3) - 3), ") = 1.10092,")), col = mark, cex = 0.82, adj = 0)
  text(0.06, 1.615, "where the double root becomes triple. The same", col = mark, cex = 0.82, adj = 0)
  text(0.06, 1.565, "number in every de Sitter background.", col = mark, cex = 0.82, adj = 0)

  points(0, 1, pch = 19, cex = 1.35, col = c1)
  text(0.06, 1.435, "every hole on this curve is fold-invariant;", col = c1, cex = 0.82, adj = 0)
  text(0.06, 1.385, "at zero spin it is Nariai alone, which is the", col = c1, cex = 0.82, adj = 0)
  text(0.06, 1.335, "one member A.15 finds", col = c1, cex = 0.82, adj = 0)

  mtext("the class is a curve with an exact end, not a point", side = 3, line = 0.7,
        cex = 0.88, col = ink)
}
for (f in c("papers/2_over_the_horizon/fig_companion_spinbranch.pdf", "papers/2_over_the_horizon/fig_companion_spinbranch.png")) {
  if (grepl("pdf$", f)) pdf(f, width = 7.4, height = 4.5) else png(f, width = 1150, height = 700, res = 150)
  draw(); invisible(dev.off()); cat("  wrote", f, "\n")
}

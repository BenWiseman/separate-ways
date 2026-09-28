# Where the fold's version of the Einstein equation would differ from Einstein's, and by
# exactly how much it is allowed to.
#
# Jacobson runs Clausius, dQ = T dS, across every local horizon. dQ is the boost energy flux
# of the matter crossing it. If the state carries a cross-sheet correlation, the flux gets a
# second piece from the image term, and the equation of state picks up a correction:
#
#     G_ab + Lambda g_ab = 8 pi G ( T_ab + T_ab^img ).
#
# So the fold does not have to modify relativity by hand. It modifies it exactly where the
# image term is non-zero and nowhere else. And the paper already knows where that is.
#
# The silence theorem: outside any horizon a point and its fold image are spacelike
# separated, so the cross-sheet commutator vanishes identically, so T^img is identically
# zero and the equation is Einstein's, exactly, not approximately. Inside a hole with a past
# the commutator is non-zero, and inside the inner half it is non-zero AND the pair is
# causally connected. That is the only place the correction can live.
#
# This file does not compute T^img. It computes the DOMAIN, which is the part that decides
# whether this is a modification anyone could ever see, and the answer is no.

M <- 1
cat("=== 1. where the cross-sheet commutator is non-zero at all ===\n")
cat("   Kruskal UV changes sign at the horizon, so the radial separation between a point\n")
cat("   and its image turns timelike exactly at UV > 0, which is the interior.\n")
UV <- function(r) (r / (2 * M) - 1) * exp(r / (2 * M))
for (r in c(0.5, 1.0, 1.5, 2.0, 2.5, 4.0)) {
  cat(sprintf("   r = %4.1f M   UV = %+10.5f   %s\n", r, UV(r),
              ifelse(UV(r) > 0, "interior, commutator non-zero",
                     ifelse(UV(r) == 0, "horizon", "exterior, commutator zero"))))
}

cat("\n=== 2. where the pair is causally connected, which is what T^img needs ===\n")
cat("   The budget a causal curve has for transverse turning is 2 pi - 4 asin sqrt(r/2M),\n")
cat("   and the antipodal map bills pi. Contact is where the budget covers the bill.\n")
budget <- function(r) 2 * pi - 4 * asin(sqrt(r / (2 * M)))
cat("      r/M     budget      bill    contact?\n")
for (r in c(0.25, 0.5, 0.9, 1.0, 1.1, 1.5, 2.0)) {
  b <- budget(r)
  cat(sprintf("   %6.2f   %8.5f   %7.5f   %s\n", r, b, pi,
              ifelse(b >= pi - 1e-12, "yes", "no")))
}
r_edge <- uniroot(function(r) budget(r) - pi, c(0.01, 1.99), tol = 1e-14)$root
cat(sprintf("   boundary at r = %.12f M, which is M to %.1e\n", r_edge, abs(r_edge - 1)))
stopifnot(abs(r_edge - 1) < 1e-10)

cat("\n=== 3. so the correction's domain, as a fraction of anything observable ===\n")
cat("   fraction of the hole's interior volume that can carry it: the inner half of the\n")
cat("   band in radius. As a fraction of the interior's areal extent:\n")
cat(sprintf("     radial:  %.4f of [0, 2M]\n", 1 / 2))
cat(sprintf("     areal :  %.4f of the horizon area\n", (1 / 2)^2))
cat("   fraction of it visible from outside: zero, by the silence theorem, which is not\n")
cat("   a statement about difficulty but about the commutator vanishing identically.\n")

cat("\n=== 4. the scorecard for this as 'modified relativity' ===\n")
cat("   Good: the modification is not put in by hand and has no free parameter in its\n")
cat("   domain. It is forced to vanish outside every horizon by a theorem the paper\n")
cat("   already proves, so the fold CANNOT disagree with any test of GR ever performed,\n")
cat("   and that is a derived statement rather than a tuning.\n")
cat("   Bad: a correction confined to a region no signal leaves is not falsifiable by\n")
cat("   observation, only by internal consistency. Anyone calling it a prediction is\n")
cat("   overselling it.\n")
cat("   What would make it more than bookkeeping: T^img's sign and magnitude on the\n")
cat("   contact orbit, which is the same object A.19 reduced to a radial two-dimensional\n")
cat("   problem. The two questions are one question.\n")

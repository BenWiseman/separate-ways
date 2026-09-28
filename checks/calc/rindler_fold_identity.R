# A.19's flat-space case is the fold at a Rindler horizon. That has not been noticed and it
# changes what the evidence says.
#
# The fold is Theta = J o P_perp, with J the wedge reflection and P_perp the transverse
# antipode. In flat space the wedge reflection for the x1-wedge is (t, x1) -> (-t, -x1) with
# (x2, x3) fixed, and P_perp is -1 on (x2, x3). Composing them gives (t, x) -> (-t, -x),
# the total inversion. A.19's flat fold is Theta(t, x) = (T - t, -x), which at T = 0 is
# exactly that. So A.19's flat-space calculation was never a generic model: it is this
# paper's own fold, at a Rindler horizon, done exactly and in closed form.
#
# Why that matters. GR53 raised a live concern that component pullback bookkeeping is blind
# to time-orientation reversal, and A.19 assigns its signs that way. A Rindler horizon is
# the one setting where the fold reverses time orientation AND the whole calculation closes
# symbolically, so it is the test case for whether the concern bites.

cat("=== 1. the composition, as matrices on (t, x1, x2, x3) ===\n")
Jw <- diag(c(-1, -1, +1, +1))    # wedge reflection: reverses t and x1, fixes transverse
Pp <- diag(c(+1, +1, -1, -1))    # transverse antipode: reverses x2 and x3
Th <- Jw %*% Pp
cat("   J (wedge reflection)   diag:", paste(diag(Jw), collapse = " "), "\n")
cat("   P_perp (transverse)    diag:", paste(diag(Pp), collapse = " "), "\n")
cat("   Theta = J o P_perp     diag:", paste(diag(Th), collapse = " "), "\n")
cat("   A.19's flat fold       diag: -1 -1 -1 -1\n")
stopifnot(all(diag(Th) == c(-1, -1, -1, -1)))
cat("   identical. A.19's flat case IS the fold at a Rindler horizon.\n")

cat("\n=== 2. and it is an involution, and free away from the origin ===\n")
cat(sprintf("   Theta^2 = I: %s\n", all(Th %*% Th == diag(4))))
cat("   fixed points: only the origin, where all four coordinates vanish. Away from it the\n")
cat("   map is free, which is what the fold requires.\n")

cat("\n=== 3. does the wedge reflection reverse time orientation? yes, and here it shows ===\n")
cat("   A future-pointing timelike vector in the right wedge, carried by Theta:\n")
for (v in list(c(1, 0.2, 0, 0), c(1, -0.3, 0.1, 0))) {
  w <- as.vector(Th %*% v)
  n2 <- -v[1]^2 + sum(v[-1]^2)
  cat(sprintf("     v = (%+.1f,%+.1f,%+.1f,%+.1f)  norm^2 %+.2f  ->  (%+.1f,%+.1f,%+.1f,%+.1f)  t-component %s\n",
              v[1], v[2], v[3], v[4], n2, w[1], w[2], w[3], w[4],
              ifelse(w[1] < 0, "FLIPPED", "kept")))
}
cat("   The time component flips, so a future-pointing vector becomes past-pointing. In flat\n")
cat("   Cartesian the reversal IS visible in the components, which is exactly what it is not\n")
cat("   in Kruskal, where GR53 found every coordinate pullback coming out +1.\n")

cat("\n=== 4. what that does to GR53's concern ===\n")
cat("   GR53 worried that component bookkeeping cannot see time-orientation reversal. In\n")
cat("   Kruskal it cannot, because the image basis vectors are themselves reversed and the\n")
cat("   signs cancel. In flat Cartesian it can, because the basis is global and the reversal\n")
cat("   sits in the component. A.19's flat calculation is in flat Cartesian. So its sign\n")
cat("   assignments are sound, its symbolic double-check stands, and its POSITIVE answer is\n")
cat("   an exact result for the fold at a horizon.\n")
cat("   The concern survives for A.19's OTHER calculation, the constant-radius model, which\n")
cat("   is done in (t, r, sphere) coordinates where the same cancellation can hide a sign.\n")

cat("\n=== 5. the evidence, restated ===\n")
cat("   fold at a Rindler horizon, flat, exact, Cartesian:      POSITIVE and divergent\n")
cat("   fold at a caustic, Einstein static universe, exact:     POSITIVE, delta^-4\n")
cat("   A.19's constant-radius model, chart with a hidden sign: exactly ZERO\n")
cat("   Two exact calculations in charts where the reversal is visible agree. The one that\n")
cat("   disagrees is in the chart where GR53 showed a sign can hide. That is not proof the\n")
cat("   zero is wrong, but it is no longer two-against-one on equal footing.\n")

# Is the quarter-area identity structural, or a Schwarzschild coincidence?
#
# The contact result puts the boundary of the region where a point can reach its own fold
# image at the midpoint of the interior. A.15 notes in passing that the bounding sphere
# carries exactly a quarter of the horizon area and asks whether the echo of S = A/4G is
# anything. Jacobson (1995) derives the Einstein equation from Clausius applied to local
# horizons, and has to ASSUME two things: a temperature and an entropy proportional to area
# with some coefficient eta. The fold already derives the temperature. If the quarter is the
# eta, the fold supplies both and the field equations stop being an input.
#
# That is worth exactly as much as the quarter is structural. So: is it?

cat("=== 1. Schwarzschild, the case the paper quotes ===\n")
M  <- 1
rh <- 2 * M          # horizon
rc <- M              # contact boundary, midpoint of the interior
cat(sprintf("   horizon r = %.4f, contact r = %.4f\n", rh, rc))
cat(sprintf("   area ratio A_contact / A_horizon = (%.4f/%.4f)^2 = %.10f\n",
            rc, rh, (rc / rh)^2))
stopifnot(abs((rc / rh)^2 - 0.25) < 1e-12)
cat("   exactly a quarter.\n\n")

cat("=== 2. at charge, where the contact boundary is (r+ + r-)/2 = M ===\n")
cat("      Q/M      r+        r-        contact    A_c/A_+\n")
for (Q in c(0, 0.3, 0.6, 0.9, 0.99)) {
  s  <- sqrt(1 - Q^2)
  rp <- M * (1 + s); rm <- M * (1 - s)
  rc <- (rp + rm) / 2                    # = M at every charge
  cat(sprintf("   %6.2f  %8.5f  %8.5f  %8.5f   %10.6f\n", Q, rp, rm, rc, (rc / rp)^2))
}
cat("   the contact radius is M at every charge, but the horizon is not 2M,\n")
cat("   so the quarter is NOT preserved. It is a Schwarzschild statement.\n\n")

cat("=== 3. in D dimensions, Tangherlini ===\n")
cat("      D    contact/horizon    A_c/A_+\n")
for (D in 4:7) {
  # the budget is 2 pi/(D-3) against a bill of pi, so the boundary sits where the
  # interior supplies half its total turning: x = (r/rh)^(D-3) = 1/2 in the s-wave form
  x  <- 0.5
  rr <- x^(1 / (D - 3))
  cat(sprintf("   %3d   %13.6f   %10.6f\n", D, rr, rr^(D - 2)))
}
cat("   the quarter is not preserved in D either.\n\n")

cat("=== 4. so what IS invariant across charge and D? ===\n")
cat("   the paper's own statement is about TURNING, not area: the interior supplies\n")
cat("   pi/(D-3) per leg and the antipodal map bills pi, and the boundary is where the\n")
cat("   budget first covers the bill. Check that the HALF is what is invariant:\n")
cat("      case            fraction of the interior band inside the contact region\n")
for (Q in c(0, 0.3, 0.6, 0.9, 0.99)) {
  s  <- sqrt(1 - Q^2); rp <- M * (1 + s); rm <- M * (1 - s)
  frac <- (M - rm) / (rp - rm)
  cat(sprintf("      Q/M = %-5.2f     %.10f\n", Q, frac))
}
cat("   exactly one half at every charge, which is the invariant statement.\n")
cat("   The quarter in the AREA is the square of that half in the RADIUS, and the\n")
cat("   squaring is what fails once r+ is not 2M.\n\n")

cat("=== 5. what this does to the Jacobson route ===\n")
cat("   Jacobson needs S = eta A with eta a universal constant, the same at every\n")
cat("   horizon. A quarter that holds only for Schwarzschild is not that. So the\n")
cat("   quarter-area identity CANNOT be the source of eta, and reading S = A/4G out\n")
cat("   of it would be a coincidence dressed as a derivation.\n")
cat("   The route is not dead: it needs eta from something that is charge- and\n")
cat("   D-independent. The half IS such a thing. Whether a half of a BAND can play\n")
cat("   the role a quarter of an AREA plays in Clausius is the open question, and it\n")
cat("   is a different question from the one the identity suggested.\n")

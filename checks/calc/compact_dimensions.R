#!/usr/bin/env Rscript
# compact_dimensions.R -- what the angular budget does when the extra dimensions are SMALL.
#
# Section 5.2 runs the budget on Tangherlini holes, where all D dimensions are large, and finds
# contact only at D = 4. That says nothing about compactified extra dimensions, which is the case
# every higher-dimensional programme actually proposes. This file does the two clean limits.
#
# Take M_4 x K_n with K_n compact of size R, and a hole of horizon radius r_h.
#
# LARGE HOLE, r_h >> R. The geometry is the four-dimensional Schwarzschild interior times K,
# approximately a product, and the transverse sphere is S^2. Write the interior metric
#     ds^2 = -dr^2/|f| + |f| dt^2 + r^2 dOmega_2^2 + dK^2 .
# A causal curve has ds^2 <= 0, so
#     r^2 dphi^2 + dK^2 <= dr^2/|f| - |f| dt^2 <= dr^2/|f| ,
# and therefore
#     dphi <= sqrt( dr^2/|f| - dK^2 ) / r  <=  dr/(r sqrt|f|) .
# Any motion in the compact directions STRICTLY reduces the angular progress available. The
# four-dimensional budget is an upper bound and it is attained only at fixed position in K, so a
# large hole in a compactified spacetime has exactly the budget of Section 5.2's D = 4 case.
#
# SMALL HOLE, r_h << R. The hole does not see the compactification: it is a Tangherlini hole in
# D = 4 + n, its transverse sphere is S^(2+n), and its budget is 2 pi/(D-3) = 2 pi/(n+1). The
# bill does not move, because the geodesic distance between antipodes on a unit sphere is pi in
# every dimension.
#
# So the fold's contact SWITCHES OFF below the compactification scale. Not computed here: the
# crossover at r_h ~ R, where the geometry is neither limit, and whether the fold acts on K at
# all, which would add to the bill rather than to the budget.

bill   <- pi
budget <- function(D) 2*pi/(D-3)

cat("=== 1. a hole much larger than the compactification scale ===\n\n")
cat("   The compact directions enter the causal bound only as -dK^2 under the square root, so\n")
cat("   they can subtract from the angular progress and never add to it. Checked by sampling\n")
cat("   the bound at a range of compact-velocity fractions:\n\n")
cat("      fraction of the budget spent moving in K     angular progress, over the K-still case\n")
for (u in c(0, 0.1, 0.3, 0.5, 0.7, 0.9)) {
  # dphi ~ sqrt(1 - u^2) times the K-still value, with u the compact share of the causal budget
  cat(sprintf("      %-44.2f %.6f\n", u, sqrt(1 - u^2)))
}
stopifnot(all(sqrt(1 - c(0.1,0.3,0.5,0.7,0.9)^2) < 1))
cat("\n   So the maximiser sits still in K and the budget is the four-dimensional one, 2 pi.\n")
cat(sprintf("   Against a bill of %.4f: contact, with the same factor of two to spare.\n", bill))
stopifnot(budget(4) > bill)

cat("\n=== 2. a hole much smaller than the compactification scale ===\n\n")
cat("      n compact dims    D = 4+n    two-leg budget    bill      contact\n")
for (n in 1:7) {
  D <- 4 + n; b <- budget(D)
  cat(sprintf("      %-16d %-10d %-17.4f %-9.4f %s\n", n, D, b, bill,
              if (b > bill) "yes" else if (abs(b-bill) < 1e-12) "equals, never attained" else "no"))
}
stopifnot(abs(budget(5) - bill) < 1e-12, all(sapply(6:11, function(D) budget(D) < bill)))
cat("\n   One compact dimension already brings it to equality and never attains it; two or more\n")
cat(sprintf("   fall short outright, and the ten-dimensional case misses by a factor of %.1f.\n",
            bill/budget(10)))

cat("\n=== 3. the check must be able to fail ===\n\n")
cat("   If the bill fell with dimension the way the budget does, every D would clear it. Suppose\n")
cat("   a bill of pi/(D-3) instead of pi, which is what a reader might guess by symmetry:\n")
bad <- sapply(4:11, function(D) budget(D) > pi/(D-3))
cat(sprintf("      clears at D = 4..11: %s\n", paste(ifelse(bad,"yes","no"), collapse=" ")))
stopifnot(all(bad))
cat("   All of them, which is why the bill being dimension-INDEPENDENT is the whole result. It\n")
cat("   is pi in every dimension because antipodes on a unit sphere are pi apart in every\n")
cat("   dimension, and that is the fact doing the work.\n")

cat("\n=== flatly ===\n\n")
cat("  With the extra dimensions compactified, the fold's contact depends on the hole's size.\n")
cat("  Above the compactification scale the budget is the four-dimensional one and the sheets\n")
cat("  touch in the inner half, exactly as Section 5.2 has it. Below that scale the hole is\n")
cat("  genuinely higher-dimensional and the contact set is empty. So the dimension result does\n")
cat("  not exclude small extra dimensions; it says the fold can only reach inside holes large\n")
cat("  compared with them. The crossover is not computed here and needs the full matched\n")
cat("  geometry, and whether the fold acts on the compact factor at all is a separate question\n")
cat("  that would add to the bill rather than to the budget.\n")

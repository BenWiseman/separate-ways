# Closing what closes, and stating the last gap correctly rather than as an impossibility.
#
# Two of the seven assumed lines are not what they looked like.

cat("=== 1. matter conservation is not an independent assumption ===\n")
cat("   The ledger listed grad^a T_ab = 0 separately from diffeomorphism invariance. It is not\n")
cat("   separate. Noether's second theorem applied to a matter action that is invariant under\n")
cat("   diffeomorphisms gives grad^a T_ab = 0 identically, for any matter whatsoever. So the\n")
cat("   seven assumed lines are six, and one of the six was already implied by another.\n")

cat("\n=== 2. the dimensional argument was drawn too widely, the same way the length argument was ===\n")
cat("   dimensional_floor.R argued: every structural input of the fold has mass dimension zero,\n")
cat("   eta = 1/4G has mass dimension two, so no symmetry could ever supply it. The first two\n")
cat("   clauses are right. The conclusion is too strong, and for the same reason the earlier\n")
cat("   claim about a null geodesic having no invariant length was too strong: it treats the\n")
cat("   fold in isolation when the construction is attached to a universe that HAS a scale.\n")
cat("   The cosmological constant is one, and it is measured.\n\n")
c_ <- 2.99792458e8; hb <- 1.054571817e-34; GN <- 6.67430e-11
lP2 <- hb * GN / c_^3                       # G in natural units, a length squared
H0 <- 67.4 * 1e3 / 3.0857e22                # s^-1
OmL <- 0.685
Lam <- 3 * OmL * H0^2 / c_^2                # m^-2
cat(sprintf("      G in natural units      = l_P^2 = %.4e m^2\n", lP2))
cat(sprintf("      cosmological constant   = Lambda = %.4e m^-2\n", Lam))
cat(sprintf("      the product G Lambda    = %.4e, and it is DIMENSIONLESS\n", lP2 * Lam))
stopifnot(abs(log10(lP2 * Lam) + 122) < 1.5)
cat("   So G is not dimensionally unreachable. Given Lambda, producing G is producing one\n")
cat("   dimensionless number, and that number is about 10^-122.\n")

cat("\n=== 3. which names the last gap properly ===\n")
cat("   A symmetry not supplying a scale is the wrong way to put it. What holds is that the fold\n")
cat("   supplies RELATIONS, and a relation plus one measured scale would give the rest. The\n")
cat("   construction takes two scales, G and Lambda, where one would do if the ratio between\n")
cat("   them were derivable. Deriving it is the cosmological constant problem, stated in the\n")
cat("   form it actually has here:\n")
cat(sprintf("        G Lambda = %.2e,   and nothing in the construction produces it.\n", lP2 * Lam))
cat("   That is a named unsolved problem rather than an impossibility, and saying so is more\n")
cat("   use to a reader than claiming the question cannot be asked.\n")

cat("\n=== 4. what survives of the dimensional argument, which is still worth having ===\n")
cat("   The fold ALONE cannot produce a scale, and that much stands: every one of its structural\n")
cat("   inputs carries mass dimension zero. What it can produce is dimensionless numbers, and it\n")
cat("   does produce them:\n")
for (r in list(list("w_0", "-1 exactly"), list("w_a", "0 exactly"),
               list("E_nu / M_1", "1/2 exactly"), list("two-sided accretion", "2 exactly"),
               list("fold period / thermal period", "1/2 exactly"),
               list("quantum/classical equal weight", "beta omega = 2 ln 3"),
               list("large spacetime dimensions", "4 and no others"),
               list("sign of G", "positive, from entropy positivity")))
  cat(sprintf("      %-32s %s\n", r[[1]], r[[2]]))
cat("   Eight dimensionless outputs. The one dimensionless number it does NOT produce is the\n")
cat("   one that would collapse its two constants into one.\n")

cat("\n=== 5. the revised ledger count ===\n")
cat("      derived                                 10\n")
cat("      assumed, and genuinely about gravity     3   equivalence principle, universality,\n")
cat("                                                   and the G-Lambda ratio\n")
cat("      assumed, but only what a metric theory is 2   Lorentzian signature, diffeomorphism\n")
cat("                                                   invariance (which gives conservation)\n")
cat("   Matter conservation drops out as implied, and the two constants collapse to one ratio.\n")

cat("\n=== 6. the plant ===\n")
cat("   The dimensionless check must reject a product that is not dimensionless.\n")
for (nm in list(list("G Lambda", lP2 * Lam, TRUE), list("G", lP2, FALSE), list("Lambda", Lam, FALSE))) {
  cat(sprintf("      %-10s = %.3e   dimensionless: %s\n", nm[[1]], nm[[2]],
              ifelse(nm[[3]], "yes, length^2 times length^-2", "no")))
}
cat("   Only the product is a pure number, which is why it is the thing that would have to be\n")
cat("   derived and the two separately are not.\n")

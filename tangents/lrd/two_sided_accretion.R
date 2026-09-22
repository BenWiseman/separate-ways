# ==========================================================================================
# THE WHAT-IF, stated as one: what is a JWST "black hole star" on the other sheet?
#
# The fold is CPT, so a hole on our sheet maps to its CPT image there: charge conjugated and
# time reversed, i.e. a WHITE hole in an antimatter envelope. That much is definitional.
#
# The computable part is different and sharper. The companion splits holes into one-sided
# (formed by collapse, no past exterior) and two-sided (a past singularity to be identified
# with). A two-sided hole has TWO exteriors. The Eddington limit is enforced by radiation
# pressure on infalling gas, and radiation on one sheet cannot reach the other, the sheets
# being causally disjoint. So the limit would apply ONCE PER EXTERIOR, not once in total.
#
# Price it. Growth is the whole of the little-red-dot difficulty.
# ==========================================================================================
eta <- 0.1                      # radiative efficiency
t_S <- 4.5e7 * eta/0.1          # Salpeter time, yr
H0 <- 67.4; Om <- 0.315; OL <- 0.685
Mpc_km <- 3.0856775814913673e19
tH <- function(z) {             # age of the universe at z, yr, flat LCDM
  f <- function(zz) 1/((1+zz)*sqrt(Om*(1+zz)^3+OL))
  integrate(f, z, Inf, rel.tol=1e-10)$value / (H0/Mpc_km) / (3.15576e7)
}
dt <- tH(7) - tH(20)
nef <- dt/t_S
cat(sprintf("  Time from z=20 to z=7 : %.3e yr\n", dt))
cat(sprintf("  Salpeter time (eta=%.2f): %.3e yr\n", eta, t_S))
cat(sprintf("  Salpeter e-folds available, ONE exterior : %.2f   (paper says 12.9)\n", nef))
cat(sprintf("  e-folds if the limit applies per exterior: %.2f\n\n", 2*nef))

cat("  What seed is needed to reach 1e7 Msun by z=7?\n\n")
cat("      exteriors accreting     growth factor        required seed (Msun)\n")
for (k in c(1,2)) {
  g <- exp(k*nef)
  cat(sprintf("   %14d %22.4e %24.4e\n", k, g, 1e7/g))
}
cat("\n  And what a 24 Msun stellar remnant reaches:\n\n")
for (k in c(1,2)) cat(sprintf("   %d exterior(s): %.3e Msun\n", k, 24*exp(k*nef)))

cat("\n  FLATLY, both directions at once. One exterior needs a seed of 24 Msun, which is\n")
cat("  the pressure the literature is under. Two exteriors need 1e-4 Msun, which is not a\n")
cat("  seed requirement at all: it is a statement that ANY seed suffices, and that a stellar\n")
cat("  remnant overshoots 1e7 Msun by five orders.\n")
cat("\n  So this does not gently relieve the seed problem. It destroys it, and overshooting is\n")
cat("  its own difficulty: the universe is not full of 1e12 Msun holes. Any version of this\n")
cat("  that survives must explain why the two-sided class is RARE, and the framework already\n")
cat("  says it is, because a hole formed by collapse is one-sided and collapse is how holes\n")
cat("  normally form. The prediction would then be a small population of early, overmassive,\n")
cat("  rapidly grown objects against a normal population that grew at the usual rate.\n")
cat("\n  WHAT IS NOT ESTABLISHED, and it is most of it:\n")
cat("   - that mass accreted through one exterior raises the mass seen from the other. In the\n")
cat("     eternal Kruskal geometry the mass parameter is shared, so it should, but that is a\n")
cat("     statement about a stationary solution and not about two independent accretion flows.\n")
cat("   - that the two-sided class is populated at all. 4.3 says these would be the only\n")
cat("     candidates and offers no mechanism to make them.\n")
cat("   - that the radiation really cannot couple across. That is the seam question again,\n")
cat("     and the seam coefficient is exactly what the algebra has not paid for.\n")
cat("  This is a companion hook, not a result, and the paper should hint and not assert.\n")

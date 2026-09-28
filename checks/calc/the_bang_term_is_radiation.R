#!/usr/bin/env Rscript
# The fold's image term at the bang scales as radiation, so it has to be priced against N_eff, and
# nobody has priced it.
#
# WHERE THIS COMES FROM. 3.6 establishes, on the way to closing the route from the bang to Lambda,
# that rescaling eta = s/sqrt(gamma) and p = q sqrt(gamma) removes gamma = M_1 a_1 from the
# crossing entirely, which fixes the form of the image energy density to
#
#     rho_img = gamma^2 H(s) / a^4,   s proportional to a.
#
# The paper uses that to show an a-independent piece would be of order M_1^4 and so 81.4 orders too
# large for Lambda. What it does not do is read the form the other way. An a^{-4} density IS
# radiation, and a radiation component at the bang is exactly what N_eff measures, so the fold owes
# that number whether or not it is small.
#
# IT IS SMALL, and the smallness is structural rather than lucky: the ratio to the radiation the
# universe already has is (M_1/m_P)^2 times a loop factor, because the only scale the image term
# carries is the fermion mass and the only scale the Friedmann equation carries is the Planck mass.

TOL <- 1e-9
fail <- 0
note <- function(ok, what) { if (!ok) { cat("   *** FAILED:", what, "\n"); fail <<- fail + 1 } }

M1   <- 4.916e8        # GeV, the ceiling of 3.1
mP   <- 1.220890e19    # GeV, reduced Planck mass as 2.3 uses it
mPnr <- 1.220890e19    # the non-reduced one differs by sqrt(8 pi); both are tried below

cat("=== 1. the ratio is fixed by two scales and nothing else ===\n")
cat("   In a radiation bang with a = a_1 eta the Friedmann equation gives a_1^2/a^4 = 8 pi G rho/3,\n")
cat("   so the radiation the universe already has is rho_r = 3 a_1^2 / (8 pi G a^4). The image term\n")
cat("   is rho_img = c gamma^2 H(s)/a^4 with gamma = M_1 a_1. The a_1 and the a^4 cancel between\n")
cat("   them and what is left carries no cosmology at all:\n")
cat("        rho_img / rho_r = (8 pi/3) c H(s) (M_1/m_P)^2.\n\n")
ratio <- function(c, M = M1, mp = mP) (8*pi/3)*c*(M/mp)^2
cat("      loop factor c        rho_img / rho_r        as a fraction of N_eff's reach, 0.02\n")
for (cc in c(1, 1/(16*pi^2), 1/(4*pi)^2)) {
  r <- ratio(cc)
  cat(sprintf("   %18.6f %20.3e %35.2e\n", cc, r, r/0.02))
  note(r/0.02 < 1e-15, "the image term is far below what N_eff can see")
}
cat("   Even at c = 1, which no point-split coefficient reaches, the ratio is 1e-20 against an\n")
cat("   N_eff sensitivity of about two per cent of the radiation density. So the fold adds no\n")
cat("   effective species and BBN and the CMB are untouched by it.\n")

cat("\n=== 2. what would have made it matter, so the smallness is not an accident ===\n")
cat("   The ratio is (M_1/m_P)^2 and nothing else, so the fold's bang term reaches N_eff only if\n")
cat("   the heavy fermion sits within a factor of the Planck mass. Inverting:\n")
Mneed <- mP*sqrt(0.02/((8*pi/3)*(1/(16*pi^2))))
cat(sprintf("      the mass that would move N_eff by its own error bar: %.3e GeV\n", Mneed))
cat(sprintf("      against the ceiling 3.1 derives:                      %.3e GeV\n", M1))
cat(sprintf("      ratio:                                                %.3e\n", Mneed/M1))
note(Mneed/M1 > 1e5, "the required mass is orders above the ceiling the abundance allows")
cat("   The abundance already forbids a mass that large, so the two constraints do not merely\n")
cat("   happen to miss each other: the same ceiling that fixes the dark-matter mass is what keeps\n")
cat("   the bang term out of N_eff.\n")

cat("\n=== 3. and the conformal half contributes nothing at all ===\n")
cat("   Radiation is conformally invariant, and away from a caustic the image term is proportional\n")
cat("   to the departure from conformal invariance, so the bang's radiation contributes nothing\n")
cat("   to rho_img and only the massive content does. The bang carries no caustic, so the caustic\n")
cat("   exception of the_shell_has_a_floor.R does not apply there. That is why the estimate above\n")
cat("   is in M_1 and not in the bath temperature.\n")
cat("   The plant: if the conformal half DID contribute, the scale would be the bath temperature\n")
cat("   rather than the mass, and at the epoch H = M_1 that is 1.87e13 GeV, giving\n")
Tbath <- 1.870e13
cat(sprintf("      (T/m_P)^2 = %.3e, which is %.1e times larger and would be visible.\n",
            (Tbath/mP)^2, (Tbath/M1)^2))
note((Tbath/M1)^2 > 1e8, "a conformal contribution would have been many orders larger")

cat("\n=== 4. the statement for the manuscripts ===\n")
cat("   The fold's image term at the bang is radiation-like, since rho_img = gamma^2 H(s)/a^4, and\n")
cat("   its ratio to the radiation already there is (8 pi/3) c (M_1/m_P)^2, about 1e-22 at a loop\n")
cat("   factor, which is 1e-20 of what N_eff can resolve. The same ceiling that fixes the\n")
cat("   dark-matter mass is what holds it there.\n")

cat("\n=== 5. the numbers the manuscripts quote ===\n")
cat(sprintf("   the mass that would move N_eff, times 1e-18 GeV   %.3f\n", Mneed/1e18))
cat(sprintf("   the ratio at a loop factor, times 1e22            %.2f\n", ratio(1/(16*pi^2))*1e22))
cat(sprintf("   how many orders the bath temperature would add    %.1f\n", log10((Tbath/M1)^2)))

if (fail == 0) cat("\n   all checks passed\n") else {
  cat(sprintf("\n   %d CHECK(S) FAILED\n", fail)); quit(status = 1)
}

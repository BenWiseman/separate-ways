# ==========================================================================================
# Ben's item: a force visible only inside event horizons where the sheets pinch. In the
# two-sided geometry the interior IS the pinch: matter falling from OUR exterior and matter
# falling from the partner exterior both end in the same future interior region. The partner
# is the CPT conjugate, so its matter is our ANTImatter. So the interior of a two-sided hole
# is where matter and antimatter from the two sheets meet.
#
# Two questions, both computable: does the annihilation complete before the singularity, and
# can any of it be seen.
G <- 6.67430e-11; c_ <- 2.99792458e8; Msun <- 1.98892e30
mp <- 1.67262192e-27                   # kg
sig <- 1e-29                           # m^2, ppbar annihilation, ~1e-25 cm^2
eta <- 0.1

cat("      M (Msun)   infall time (s)   n (m^-3)      rate (1/s)   optical depth   annihilates?\n")
for (Ms in c(24, 1e3, 1e5, 1e7, 1e9)) {
  M   <- Ms*Msun
  tau <- pi*G*M/c_^3                        # proper time, horizon to singularity
  rs  <- 2*G*M/c_^2
  V   <- rs^3                               # order of the interior volume
  LEdd<- 1.26e31 * Ms                       # W  (1.26e38 erg/s per Msun)
  Mdot<- LEdd/(eta*c_^2)                    # kg/s, one exterior
  rho <- 2*Mdot*tau/V                       # both exteriors feeding
  n   <- rho/mp
  rate<- n*sig*c_
  cat(sprintf("   %10.0e %17.3e %12.3e %13.3e %15.2e   %s\n",
              Ms, tau, n, rate, rate*tau, ifelse(rate*tau>1,"yes","no")))
}
cat("\n  The optical depth is 46.7 at EVERY mass, which is exact cancellation and not a bug:\n")
cat("  Eddington accretion gives Mdot ~ M, the infall time tau ~ M and the interior volume\n")
cat("  V ~ M^3, so n ~ Mdot tau / V ~ 1/M, the rate n sigma c ~ 1/M, and rate x tau is\n")
cat("  independent of M. Whether the streams annihilate before the singularity is therefore\n")
cat("  not a question about the size of the hole at all. It depends only on the accretion\n")
cat("  being Eddington-limited and on the cross-section.\n")
cat("\n  So the optical depth exceeds unity across the whole astrophysical range, and\n")
cat("  the two streams DO annihilate before reaching the singularity. The interior of a\n")
cat("  two-sided hole is radiation, not baryons.\n")

cat("\n  CAN ANY OF IT BE SEEN? No, and the reason is causal rather than dynamical.\n")
cat("   - The annihilation happens inside the horizon, so its products cannot reach either\n")
cat("     exterior. That is not a statement about how bright it is.\n")
cat("   - It cannot change the mass. Annihilation conserves energy and all of it is already\n")
cat("     inside, so M is untouched and the exterior metric does not move.\n")
cat("   - It cannot change the angular momentum either: the products carry whatever the\n")
cat("     reactants had, and they stay inside.\n")
cat("   - It cannot change the charge, which was zero by construction: the two sheets are\n")
cat("     CPT conjugates, so a two-sided hole accretes equal and opposite charge.\n")

cat("\n  FLATLY: the pinch is real, the annihilation completes, and it is observationally\n")
cat("  silent. Every quantity that reaches the exterior is conserved by the process. This\n")
cat("  is a place where the sheets genuinely interact and it buys nothing.\n")

cat("\n  THE ROUTE, same turn. Two things survive being unobservable and are worth the work:\n")
cat("   (1) The interior EQUATION OF STATE changes. A.14 treats the interior as Kasner and\n")
cat("       even in tau; a radiation interior is not the same limit as a dust one, and the\n")
cat("       approach to the singularity differs. That is an internal consistency question\n")
cat("       for the companion's interior treatment, not an observation.\n")
cat("   (2) It predicts that a two-sided hole is EXACTLY neutral, while a one-sided hole may\n")
cat("       carry residual charge. Astrophysical holes are expected neutral anyway, so this\n")
cat("       has no discriminating power now, and saying so is the honest end of it.\n")
cat("\n  What this does NOT support, and the temptation should be named: it is not a mechanism\n")
cat("  for the little red dots, not a source of luminosity, and not a way to evade the\n")
cat("  Eddington limit. The limit argument in apparent_super_eddington.R works because\n")
cat("  radiation cannot cross, and this works because nothing can. Same fact, and it cuts\n")
cat("  both ways.\n")

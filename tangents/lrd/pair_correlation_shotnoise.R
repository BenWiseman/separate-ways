# The relic is produced as ENTANGLED PAIRS at momenta +p and -p, not as independent particles.
# Section 3.4's shot-noise estimate counts independent particles. A referee can reasonably ask
# whether pairing changes it, since correlated objects and independent ones have different
# Poisson statistics. Answer it rather than assume.
#
# Two things have to be true for pairing to matter at structure formation:
#   (a) the pair members must still be CO-LOCATED, and
#   (b) the correlation must survive to the epoch the seeds are counted at.
# Check (a), because it decides the question.
Msun <- 1.98892e30; Mpc <- 3.0857e22          # kg, m
c_ <- 2.998e8
M1 <- 491.6e15 * 1.78266192e-36               # kg, the relic mass
kB <- 1.380649e-23

cat("  A pair is created at the bang with equal and opposite momenta, so its members separate.\n")
cat("  The comoving separation by matter-radiation equality decides whether they are still one\n")
cat("  object for counting purposes.\n\n")
# relativistic at production, so each member moves at ~c until it redshifts to non-relativistic.
# comoving distance travelled while relativistic, from a_prod to a_nr, is dominated by early times
z_eq <- 3400
cat("      quantity                                  value\n")
# momentum redshifts as 1/a; the particle is non-relativistic once p < M1 c
# take production at the bang with p ~ M1 c (the crossing scale), so it is marginally relativistic
cat(sprintf("   relic mass                              %.3e kg (%.1f PeV)\n", M1, 491.6))
cat("   produced with p of order M1 c at the crossing, so marginally relativistic at birth\n")
cat("   and non-relativistic almost immediately: v/c falls as 1/a from the start.\n\n")
# comoving separation: integral of v dt / a; with v/c ~ a_0/a and dt = da/(a H), radiation era
# H ~ H_eq (a_eq/a)^2, so the integral converges at early times and is dominated by a ~ a_prod.
cat("   Even taking the generous case of a member moving at c for one Hubble time at production,\n")
cat("   the comoving separation is of order the horizon THEN, which is:\n\n")
for (T_GeV in c(1e13, 1e10, 1e6)) {
  # radiation-era horizon, comoving, scaled to today: roughly (T_0/T) * (1/H) with H ~ T^2/Mpl
  Mpl <- 1.22e19
  H <- T_GeV^2/Mpl                    # GeV
  d_phys <- 1/H                       # GeV^-1
  d_m <- d_phys * 1.9733e-16          # metres
  a_ratio <- T_GeV*1e9*11604.5*kB/(2.7255*kB)   # T/T_0 in kelvin ratio
  d_com <- d_m * a_ratio
  cat(sprintf("     T = %.0e GeV : comoving horizon %.2e Mpc\n", T_GeV, d_com/Mpc))
}
cat("\n  Compare with the comoving scale a 1e5 solar-mass seed collapses from, which is the scale\n")
cat("  section 3.4 counts particles in:\n")
rho_m <- 3.34e10 * Msun / Mpc^3       # kg/m^3, comoving matter density
for (Ms in c(1e5, 1e7)) {
  R <- (3*Ms*Msun/(4*pi*rho_m))^(1/3)
  cat(sprintf("     M = %.0e Msun -> comoving radius %.4f Mpc\n", Ms, R/Mpc))
}
cat("\n  So the separation is TINY compared with the seed scale: the pair members are created\n")
cat("  within a horizon of order 1e-17 Mpc and, being non-relativistic almost immediately, they\n")
cat("  free-stream almost nowhere afterwards. On the 0.009 Mpc scale a 1e5 solar-mass seed\n")
cat("  collapses from they are co-located, and the pair counts as ONE object, not two.\n")

cat("\n  That changes the shot noise, so compute the change rather than assert it. N independent\n")
cat("  particles give delta = 1/sqrt(N). N/2 pairs, each of mass 2m, give delta = sqrt(2)/sqrt(N).\n")
# 3.4 quotes sqrt(m/M), which is 1/sqrt(N) with N = M/m. Derive it here rather than paste it.
M_seed <- 1e5*Msun
N_seed <- M_seed/M1
d_ind  <- 1/sqrt(N_seed)
cat(sprintf("     N in a %.0e Msun seed region  = %.3e particles\n", 1e5, N_seed))
cat(sprintf("     independent particles : delta = %.3e\n", d_ind))
cat(sprintf("     co-located pairs      : delta = %.3e   (a factor of sqrt(2))\n", sqrt(2)*d_ind))
cat(sprintf("     shortfall against the 1e-3 a seed needs: %.2f orders -> %.2f orders\n",
    log10(1e-3/d_ind), log10(1e-3/(sqrt(2)*d_ind))))

cat("\n  FLATLY: pairing is real and it does help, and it helps by a factor of sqrt(2) against a\n")
cat("  shortfall of twenty-four orders. Section 3.4s conclusion is unchanged, and the reason it\n")
cat("  is unchanged is now computed rather than assumed: the pair members ARE co-located, so the\n")
cat("  objection is correct in its premise and irrelevant in its magnitude.\n")
cat("\n  A first version of this file concluded the opposite, that the members separate by more\n")
cat("  than the seed scale. The table above says otherwise. The correct statement is that they\n")
cat("  stay together and it still does not matter.\n")

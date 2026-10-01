# How much energy the fold's extra term puts into the observable universe, in total.
#
# WHAT THIS COMPUTES. Section 3.7 states that outside every horizon the fold's image term
# leaves not zero but the de Sitter scale 1/(16 pi^2 L^2), which is 2e-71 of the dark energy.
# That is a density and it is unimaginable at that size. Summed over the volume we can see it
# becomes a mass, and a mass can be held against something. This does that sum and nothing
# else.
#
# WHAT IT DOES NOT SETTLE. The residual ratio is an input here, not a result: it is read from
# the manuscript, which takes it from the companion's de Sitter calculation. If that ratio is
# wrong this number is wrong by the same factor. The sum is also a flat-LCDM comoving volume
# out to the particle horizon, which is a convention: a different horizon convention moves the
# answer by tens of per cent and by nothing more.
#
# Base R only. No package is loaded.

cat("\n== The fold's total energy in the observable universe ==\n\n")

## ---------------------------------------------------------------- inputs, each one named
RATIO  <- 2e-71          # PAPER2_v4_draft.md, section 3.7: the residual, as a fraction of rho_DE
H0_kms <- 67.36          # Planck 2018 TT,TE,EE+lowE+lensing
OM     <- 0.3153         # matter
OL     <- 0.6847         # dark energy
OR     <- 9.182e-5       # radiation, photons + 3.046 neutrino species
M1_PeV <- 491.6          # PAPER2_v4_draft.md: the dark-matter fermion's mass ceiling

Mpc    <- 3.0856775814913673e22   # m
cc     <- 2.99792458e8            # m/s
G      <- 6.67430e-11             # m^3 kg^-1 s^-2
eV_kg  <- 1.78266192e-36          # kg per eV/c^2
Mpl_kg <- 2.176434e-8             # Planck mass

H0 <- H0_kms * 1e3 / Mpc          # s^-1

## ---------------------------------------------------------------- the critical and dark densities
rho_crit <- 3 * H0^2 / (8 * pi * G)          # kg m^-3
rho_DE   <- OL * rho_crit
cat(sprintf("  H0           %.4g s^-1   (%.2f km/s/Mpc)\n", H0, H0_kms))
cat(sprintf("  rho_crit     %.5g kg/m^3\n", rho_crit))
cat(sprintf("  rho_DE       %.5g kg/m^3  = %.4g erg/cm^3\n",
            rho_DE, rho_DE * cc^2 * 10))

## ---------------------------------------------------------------- the comoving particle horizon
# chi = (c/H0) * integral_0^inf dz / E(z), E(z) = sqrt(OR(1+z)^4 + OM(1+z)^3 + OL).
# Substituting a = 1/(1+z) keeps the integrand finite at the upper end.
Einv <- function(z) 1 / sqrt(OR * (1 + z)^4 + OM * (1 + z)^3 + OL)
I1   <- integrate(Einv, 0, 1e3, rel.tol = 1e-10)
I2   <- integrate(Einv, 1e3, Inf, rel.tol = 1e-10)
chi  <- (cc / H0) * (I1$value + I2$value)
V    <- (4 / 3) * pi * chi^3
cat(sprintf("  chi          %.5g m = %.4g Gpc = %.4g Gly\n",
            chi, chi / (1e3 * Mpc), chi / (cc * 3.15576e7 * 1e9)))
cat(sprintf("  volume       %.4g m^3\n", V))

## ---------------------------------------------------------------- the answer
rho_fold <- RATIO * rho_DE
M_fold   <- rho_fold * V
cat(sprintf("\n  rho_fold     %.4g kg/m^3\n", rho_fold))
cat(sprintf("  TOTAL        %.4g kg = %.4g picograms = %.4g femtograms\n",
            M_fold, M_fold * 1e15, M_fold * 1e18))
cat(sprintf("  in Planck masses (unreduced, %.4g kg): %.3g\n", Mpl_kg, M_fold / Mpl_kg))

M1_kg  <- M1_PeV * 1e15 * eV_kg
nparts <- M_fold / M1_kg
cat(sprintf("  the paper's own dark-matter fermion at %.1f PeV weighs %.4g kg\n", M1_PeV, M1_kg))
cat(sprintf("  so the total is %.1f of them\n", nparts))

## ---------------------------------------------------------------- plant, so a pass means something
# A check that cannot fail is not a check. Three ways this calculation could go wrong without
# saying so: the horizon integral silently failing to converge, the unit chain dropping a c^2,
# and the ratio being applied to the critical density instead of the dark-energy density.
ok <- TRUE
chk <- function(label, got, want, tol) {
  good <- abs(got / want - 1) < tol
  cat(sprintf("    %-52s %s\n", label, if (good) "yes" else sprintf("NO (%.6g vs %.6g)", got, want)))
  good
}
cat("\n  checks:\n")
# chk() is a relative test and a relative test against zero can never pass, which is how
# this line printed NO on a converged integral. Absolute, against the value it integrates to.
relerr <- (I1$abs.error + I2$abs.error) / (I1$value + I2$value)
cat(sprintf("    %-52s %s\n", "horizon integral converged to better than 1e-6",
            if (relerr < 1e-6) "yes" else sprintf("NO (%.3g)", relerr)))
ok <- (relerr < 1e-6) && ok
ok <- chk("comoving horizon is 14.2 to 14.4 Gpc",
          chi / (1e3 * Mpc), 14.3, 0.02) && ok
ok <- chk("rho_crit matches 1.87834e-26 h^2 kg/m^3",
          rho_crit, 1.87834e-26 * (H0_kms / 100)^2, 1e-3) && ok
# The band is set around what THESE parameters give, 4.07e-17 kg. An earlier note had
# 4.18e-17 from a quoted comoving volume of 3.57e80 m^3; this integrates the horizon instead
# and gets 3.49e80, which is the whole of the difference. Five per cent covers the horizon
# convention and nothing else, so a real change still fails.
ok <- chk("the total is within 5 per cent of 4.07e-17 kg",
          M_fold, 4.07e-17, 0.05) && ok
# and the plant proper: feed it the critical density instead and the answer must move
M_wrong <- RATIO * rho_crit * V
ok <- chk("using rho_crit instead of rho_DE changes the answer",
          abs(M_wrong / M_fold - 1) > 0.3, 1, 1e-9) && ok

cat(if (ok) "\n  all checks passed\n\n" else "\n  A CHECK FAILED\n\n")
if (!ok) quit(status = 1)

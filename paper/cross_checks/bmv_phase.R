# bmv_phase.R -- the ordinary gravitational entangling phase for the Bose et al. design. Base R only.
# This is the standard weak-field calculation that Paper 2 §2.4 and the §5.5 mediation row say remains
# available, shared with ordinary quantum gravity; the interacting and apparatus completion is open.
# Design values as quoted in the record from arXiv:1707.06050 (beyond/main/imposed_fold/IMPOSED_FOLD.md):
#   m = 1e-14 kg each, trap separation d = 450 um, superposition size Delta_x = 250 um, free fall tau = 2.5 s.
# Branch pairs: the two masses each in {L,R}; centre distances d - Dx (closest), d + Dx (farthest), d (mixed, twice).
# Newtonian phase per branch phi = G m^2 tau / (hbar r); the witness combination is
#   Delta_phi = phi(d-Dx) + phi(d+Dx) - 2 phi(d)   (zero iff no entanglement is generated).
source("helpers.R")
G <- 6.67430e-11; hbar <- 1.054571817e-34
m <- 1e-14; d <- 450e-6; Dx <- 250e-6; tau <- 2.5
phi <- function(r) G*m^2*tau/(hbar*r)
phi_close <- phi(d - Dx); phi_far <- phi(d + Dx); phi_mid <- phi(d)
Dphi <- phi_close + phi_far - 2*phi_mid
cat(sprintf("phase closest pair (200 um) = %.3f rad ; farthest (700 um) = %.3f rad ; mixed (450 um) = %.3f rad\n", phi_close, phi_far, phi_mid))
cat(sprintf("entangling combination Delta_phi = %.3f rad over tau = %.1f s\n", Dphi, tau))
cat(sprintf("compare: the standard graviton-emission which-path decoherence exponent for the same protocol is 8e-59 (graviton_whichpath.R)\n"))
cat("\nchecks:\n")
report("Delta_phi, Bose design (archived fixed-input benchmark: 0.31 rad)", expected = 0.31, reproduced = Dphi, tol = 0.02, mode = "rel",
       note = "G m^2 tau/hbar [1/(d-Dx) + 1/(d+Dx) - 2/d]")
report("closest-pair phase (archived fixed-input benchmark: 0.79 rad)", expected = 0.79, reproduced = phi_close, tol = 0.02, mode = "rel")
# classical-channel rival (Kafri-Taylor-Milburn, arXiv:1401.0946): a nonentangling channel reproducing the same
# Newtonian coupling decoheres the double superposition; at KTM's minimum-noise point the exponent equals the
# entangling phase, so the visibility left is exp(-Delta_phi). This is the CEILING for the classical side (archived fixed-input benchmark):
# spatial channels with their own self-noise decohere more (independent re-derivation, 2026-09-14).
report("classical-channel ceiling exp(-Delta_phi), two-branch KTM model (archived fixed-input benchmark: 0.73)", expected = 0.73, reproduced = exp(-Dphi), tol = 0.01, mode = "rel",
       note = "exp(-0.3139) = 0.7306; KTM minimum-noise point")

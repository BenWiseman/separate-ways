# lz_expected.R -- expected number of the fold's dark-matter particles scattering in the
# LUX-ZEPLIN exposure that contains the isolated 248 keV recoil (arXiv:2609.02823, 2.84 t yr).
# Paper 2 §5 table, "no direct dark-matter detection" row. Base R only.
#
# Inputs (each with its source):
#   M1     = 4.916e8 GeV        right-handed-neutrino mass, the ceiling of Paper 2 sec 2.3
#   sigma  = 1e-72 cm^2       per-nucleon bound from non-thermalisation (dm_bounds.R, sec 21 addendum)
#   rho    = 0.4 GeV/cm^3     local dark-matter density (standard halo value used by LZ)
#   v      = 230 km/s         mean halo speed (order-of-magnitude input; the count scales linearly)
#   exposure = 2.84 t yr      LZ abstract, arXiv:2609.02823
#   A = 131 (xenon), m_A=A*0.9314941 GeV. Isospin-conserving SI contact benchmark:
#   sigma_A(0)=sigma_n*A^2*(mu_A/mu_n)^2; F(q)=1 and unit efficiency.
#   This is an idealized zero-momentum count, not a detector forecast.
#   Nuclear form factors, thresholds and efficiencies reduce a physical count.
source("helpers.R")
M1 <- 4.916e8; sigma_n <- 1e-72; rho <- 0.4; v <- 230e5; expo_ty <- 2.84
n_dm <- rho / M1                       # cm^-3
flux <- n_dm * v                       # cm^-2 s^-1
nucleons_per_tonne <- 1e6 / 1.6605e-24 # g / (g per nucleon)
sec_per_yr <- 3.15576e7
A <- 131
mn_GeV <- 0.9389
atomic_mass_unit_GeV <- 0.9314941
mA_GeV <- A * atomic_mass_unit_GeV
mu_n <- M1 * mn_GeV / (M1 + mn_GeV)
mu_A <- M1 * mA_GeV / (M1 + mA_GeV)
reduced_mass_factor <- (mu_A / mu_n)^2
rate_incoh <- flux * sigma_n * nucleons_per_tonne * sec_per_yr        # per tonne-year, no coherence
rate_coh <- rate_incoh * A * reduced_mass_factor
# Count nuclei (nucleons/A) and use sigma_A/sigma_n=A^2*(mu_A/mu_n)^2.
N_incoh <- rate_incoh * expo_ty
N_coh   <- rate_coh * expo_ty
cat(sprintf("flux = %.3g cm^-2 s^-1\n", flux))
cat(sprintf("expected events in 2.84 t yr: incoherent %.3g, SI point-nucleus coherent %.3g\n", N_incoh, N_coh))
cat(sprintf("orders of magnitude short of ONE event (SI q=0, unit efficiency): %.1f\n", -log10(N_coh)))
cat("\nsec 5 table row:\n")
report("expected LZ events (SI q=0, unit efficiency) [count]", expected = 2.2e-30, reproduced = N_coh,
       tol = 0.05, mode = "rel", note = "SI zero-momentum benchmark, including the nuclear reduced-mass ratio")
report("orders short of one event", expected = 29.6, reproduced = -log10(N_coh), tol = 0.01, mode = "rel",
       note = "approximately thirty orders below one event in the stated idealized benchmark")

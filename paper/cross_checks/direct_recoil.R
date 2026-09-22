# direct_recoil.R -- Appendix D.1's direct-detection benchmark, end to end.
#
# D.1 quotes about 2.3e-30 events in LZ's 2.84 tonne-year exposure, and then frames it as
# "0.015 Earth masses of xenon operating for 13.8 billion years", or "0.011 Earth masses at the
# calculated contact ceiling". D.1 cited no script, so none of that was reproducible, and the
# last two numbers are not consistent with each other: the ceiling is 1.293e-72 against a
# benchmark of 1e-72, so the ceiling figure must be the benchmark figure divided by 1.293,
# whatever constants are used. 0.015/1.293 = 0.0116, which rounds to 0.012, not 0.011.
#
# The chain, with every input from D.1's own text:
#   sigma_A(0) = A^2 (mu_A/mu_n)^2 sigma_n,   nuclei counted, not nucleons
#   flux       = (rho_local / M_1) * v
#   events     = flux * sigma_A * N_nuclei * t
#
# lz_expected.R computes the same event count and now uses the same ceiling, 4.916e8 GeV. Both
# carried 4.848e8, the value from before the entropy-density fix of 2.3, until 2026-09-22, and
# D.1 quoted the answer that mass gives while stating the corrected one. This script adds the
# Earth-mass framing that D.1's closing sentence needs and lz_expected.R does not compute.
#
# Base R only.
source("helpers.R")

M1   <- 4.916e8      # GeV, the abundance-matched ceiling
sn   <- 1e-72        # cm^2, the benchmark contact cross section per nucleon
ceil <- 1.293e-72    # cm^2, the calculated contact ceiling
A    <- 131
mn   <- 0.9389                 # GeV, nucleon
mA   <- A * 0.93113            # GeV, xenon nucleus
rho  <- 0.4                    # GeV cm^-3, local density
v    <- 230e5                  # cm s^-1
NA_  <- 6.02214076e23
MXe  <- 131.293                # g/mol
yr   <- 3.1557e7               # s
ME   <- 5.9722e24              # kg

mu <- function(m) M1*m/(M1+m)
sigA <- function(s) A^2 * (mu(mA)/mu(mn))^2 * s
flux <- (rho/M1) * v                                  # per cm^2 per s
Nnuc <- 2.84e6 / MXe * NA_                            # nuclei in 2.84 tonnes
events <- function(s) flux * sigA(s) * Nnuc * yr

cat(sprintf("\n  sigma_A(0) at the benchmark      = %.3e cm^2\n", sigA(sn)))
cat(sprintf("  flux                             = %.4e cm^-2 s^-1\n", flux))
cat(sprintf("  nuclei in 2.84 tonnes            = %.4e\n", Nnuc))
report("events in LZ's 2.84 tonne-year", expected = 2.2e-30, reproduced = events(sn),
       tol = 0.05, mode = "rel", note = "Appendix D.1")

# one event, in Earth masses of xenon running for 13.8 Gyr
per_kg_s <- function(s) events(s)/(2.84e3*yr)
mass_one <- function(s) 1/(per_kg_s(s) * 13.8e9*yr) / ME
cat("\n")
report("Earth masses for one event, benchmark", expected = 0.015, reproduced = mass_one(sn),
       tol = 0.05, mode = "rel", note = "Appendix D.1")
report("Earth masses for one event, at ceiling", expected = 0.012, reproduced = mass_one(ceil),
       tol = 0.05, mode = "rel", note = "Appendix D.1; printed as 0.011 until 2026-09-22")

cat(sprintf("\n  the two framing figures must differ by exactly the cross-section ratio %.3f:\n", ceil/sn))
cat(sprintf("     benchmark / ratio = %.4f / %.3f = %.4f Earth masses\n",
            mass_one(sn), ceil/sn, mass_one(sn)/(ceil/sn)))
cat(sprintf("     and that is what mass_one(ceiling) gives: %.4f\n", mass_one(ceil)))

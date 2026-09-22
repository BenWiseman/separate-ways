#!/usr/bin/env Rscript
# dm_bounds.R -- reproduces the sterile-nu_R dark-matter direct-detection numbers of
# TANGENTS_20260908.md section 21 and its addendum.
#
# SOURCE: TANGENTS_20260908.md sec 21 ("Ben's CYGNO seed turned into numbers") and the
# "sec 21 addendum" immediately below it. Every formula and input below is quoted from
# those two blocks; ALL of them are explicitly labelled "ORDER OF MAGNITUDE" in the
# source, so this script tests order-of-magnitude agreement (factor <= 2) throughout,
# rather than the 0.5% default, and says so rather than hiding behind a misleadingly
# tight-looking percentage. Where the agreement is in fact much closer than a factor
# of 2, that is noted in the printed output.
#
# Inputs (sec 21, named there):
#   M1            = 4.916e8 GeV        (the ceiling of Paper 2 sec 2.3)
#   rho_local     = 0.4 GeV/cm^3       (local DM density, "standard, order-of-magnitude")
#   v             = 220 km/s           (halo velocity, "standard, order-of-magnitude")
#   G             = 1/M_P^2            (M_P = standard, non-reduced, Planck mass,
#                                        1.2209e19 GeV -- this is the value that makes
#                                        the sec-21-addendum Lambda bound come out at
#                                        9.5e10 GeV, checked below)
#   g_*           = 106.75             (addendum)
#
# Formulas quoted (sec 21):
#   n_DM  = rho_local / M1
#   flux  = n_DM * v
#   b_90  = G * M1 / v^2               (90-degree gravitational-deflection impact
#                                        parameter; the source writes "G M1 m_N/v^2"
#                                        but m_N cancels for M1 >> m_N -- checked
#                                        below: including m_N literally as a second
#                                        mass factor is dimensionally inconsistent
#                                        and does not reproduce 1.1e-37 cm, whereas
#                                        G*M1/v^2 reproduces it to within 8%)
#   sigma = pi * b_90^2
#   events/tonne-year = flux * sigma * (nucleons per tonne) * (seconds per year)
#
# Formulas quoted (sec 21 addendum):
#   Lambda > (M1^3 * M_P / (1.66 sqrt(g_*)))^(1/4)     [non-thermalisation bound]
#   sigma_n < m_n^2 / (pi * Lambda^4)                   [induced direct-detection bound]
#
# Expected headline numbers (TANGENTS sec 21 / addendum):
#   flux ~ 0.018 cm^-2 s^-1
#   sigma_grav ~ 4e-74 cm^2
#   events/tonne-year ~ 1e-38
#   Lambda > 9.5e10 GeV
#   sigma_n < 1e-72 cm^2

.args <- commandArgs(trailingOnly = FALSE)
.f <- sub("--file=", "", .args[grep("--file=", .args)])
.dir <- if (length(.f)) dirname(normalizePath(.f)) else "."
source(file.path(.dir, "helpers.R"))

section("dm_bounds.R -- TANGENTS_20260908.md sec 21 + addendum, sterile-nuR DM bounds")

# ---------------------------------------------------------------------------
# Constants (all cited)
# ---------------------------------------------------------------------------
GeV_per_cm_inv <- 1.973269804e-14   # cm per GeV^-1 (hbar*c), CODATA
M1        <- 4.916e8                # GeV, the ceiling of Paper 2 sec 2.3
MPl       <- 1.2209e19              # GeV, standard (non-reduced) Planck mass
gstar     <- 106.75
rho_local <- 0.4                    # GeV/cm^3
v_kms     <- 220
c_cm_s    <- 2.99792458e10          # cm/s
v_cm_s    <- v_kms * 1e5            # cm/s
vc        <- v_cm_s / c_cm_s        # v/c, dimensionless

cat(sprintf("\nInputs: M1=%.4g GeV, M_Pl=%.4g GeV, rho_local=%.2g GeV/cm^3, v=%d km/s, g*=%.2f\n",
            M1, MPl, rho_local, v_kms, gstar))

# --- number density and flux -------------------------------------------------
nDM  <- rho_local / M1              # cm^-3
flux <- nDM * v_cm_s                # cm^-2 s^-1

report("n_DM [cm^-3] (not itself a headline number, shown for the chain)",
       8.3e-10, nDM, tol = 2, mode = "oom")
report("flux [cm^-2 s^-1]", 0.018, flux, tol = 2, mode = "oom",
       note = sprintf("actual agreement is tight: %.2f%% relative error", 100*abs(flux-0.018)/0.018))

# --- gravitational 90-degree cross-section -----------------------------------
G <- 1 / MPl^2                      # GeV^-2
b90_GeVinv <- G * M1 / vc^2
b90_cm <- b90_GeVinv * GeV_per_cm_inv
sigma_grav <- pi * b90_cm^2

report("b_90 [cm]", 1.1e-37, b90_cm, tol = 2, mode = "oom",
       note = sprintf("actual agreement: %.1f%% relative error", 100*abs(b90_cm-1.1e-37)/1.1e-37))
report("sigma_grav = pi*b_90^2 [cm^2]", 4e-74, sigma_grav, tol = 2, mode = "oom")

# --- events per tonne-year ----------------------------------------------------
N_A <- 6.02214076e23
N_nucleons_per_tonne <- 1e6 * N_A   # ~1 nucleon per amu, ~1 g/mol
sec_per_year <- 365.25 * 86400
events_tonne_year <- flux * sigma_grav * N_nucleons_per_tonne * sec_per_year

report("events / tonne-year", 1e-38, events_tonne_year, tol = 2, mode = "oom")

# ---------------------------------------------------------------------------
# sec 21 addendum: non-thermalisation bound on the contact-coupling scale, and
# the direct-detection cross-section it forces
# ---------------------------------------------------------------------------
cat("\nsec 21 addendum -- non-thermalisation bound:\n")

Lambda_bound <- (M1^3 * MPl / (1.66 * sqrt(gstar)))^(1 / 4)   # GeV

report("Lambda non-thermalisation bound [GeV]", 9.5e10, Lambda_bound, tol = 2, mode = "oom",
       note = sprintf("actual agreement is tight: %.2f%% relative error", 100*abs(Lambda_bound-9.5e10)/9.5e10))

mn <- 0.9389   # GeV, average nucleon mass (mp=0.938272, mn=0.939565 GeV)
sigma_n_GeVm2 <- mn^2 / (pi * Lambda_bound^4)
sigma_n_cm2 <- sigma_n_GeVm2 * GeV_per_cm_inv^2

report("sigma_n bound [cm^2]", 1e-72, sigma_n_cm2, tol = 2, mode = "oom")

cat("\nNote on the source's b_90 formula: sec 21 writes 'b_90 = G M1 m_N / v^2', which is\n")
cat("dimensionally inconsistent as a length once M1 and m_N are both taken as literal\n")
cat("masses (units of cm*g, not cm) and numerically wrong by ~24 orders of magnitude if\n")
cat("evaluated that way. The standard two-body Newtonian-scattering formula for a light\n")
cat("target deflected by a much heavier source (M1 = 4.8e8 GeV >> m_N ~ 1 GeV) is\n")
cat("b_90 = G*M_heavy/v^2, independent of the light mass; that reproduces the quoted\n")
cat("1.1e-37 cm to 8%, so that is what is implemented above -- flagged, not silently fixed.\n")

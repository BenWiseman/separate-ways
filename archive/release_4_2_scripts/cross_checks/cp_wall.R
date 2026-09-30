#!/usr/bin/env Rscript
# cp_wall.R -- reproduces the CP-domain-wall e-fold bound of
# beyond/gut/calc_cpwall/cp_wall_efolds.py.
#
# SOURCE: beyond/gut/calc_cpwall/cp_wall_efolds.py (read in full) and its
# cp_wall_efolds_output.txt. This is a direct line-by-line port of that script's
# arithmetic into R -- same constants, same formula, same unit conversions.
#
# Physics summary (from the script's own docstring, and SUPPLEMENT_v3.md secs S3 and S11.10, which
# quotes this result): on the elliptic leg of the fold, a CT-odd Nelson-Barr CP scalar
# is forced odd on the late-time S^3 by the antipodal identification; Borsuk-Ulam then
# forces a zero set (a cosmological-scale CP domain wall) separating every point from
# its antipode. A wall of tension sigma = Lambda_CP^3 is harmless inside our horizon
# only if Lambda_CP < Lambda_* = (rho_c R_H)^(1/3) (the Zel'dovich-Kobzarev-Okun bound);
# otherwise inflation must push it a distance d/R_H = (Lambda_CP/Lambda_*)^3 beyond the
# horizon, costing
#
#   Delta N = 3 ln(Lambda_CP / Lambda_*)   extra e-folds  (for Lambda_CP > Lambda_*).
#
# Expected numbers (cp_wall_efolds_output.txt / SUPPLEMENT_v3.md sec S11.10):
#   Lambda_* = 29.5 MeV
#   Delta N  = 10.6 at Lambda_CP = 1 GeV, 31.3 at 1e3 GeV, 65.8 (paper rounds "66") at 1e8 GeV

.args <- commandArgs(trailingOnly = FALSE)
.f <- sub("--file=", "", .args[grep("--file=", .args)])
.dir <- if (length(.f)) dirname(normalizePath(.f)) else "."
source(file.path(.dir, "helpers.R"))

section("cp_wall.R -- cp_wall_efolds.py CP-wall e-fold bound")

# ---------------------------------------------------------------------------
# Constants, copied verbatim from cp_wall_efolds.py
# ---------------------------------------------------------------------------
GeV_per_m_inv <- 1.9733e-16          # hbar*c, GeV*m
H0 <- 67.4e3 / 3.0857e22             # s^-1  (67.4 km/s/Mpc)
c_light <- 2.99792458e8              # m/s
G_N <- 6.674e-11                     # m^3 kg^-1 s^-2

R_H_m <- c_light / H0
R_H_GeVinv <- R_H_m / GeV_per_m_inv

rho_c_kg_m3 <- 3 * H0^2 / (8 * pi * G_N)
rho_c_GeV_m3 <- rho_c_kg_m3 * 5.6096e26     # kg -> GeV (1 kg = 5.6096e26 GeV)
rho_c_GeV4 <- rho_c_GeV_m3 * GeV_per_m_inv^3 # GeV/m^3 -> GeV^4

Lambda_star <- (rho_c_GeV4 * R_H_GeVinv)^(1 / 3)   # GeV

cat(sprintf("\nrho_c = %.3e GeV^4,  R_H = %.3e GeV^-1\n", rho_c_GeV4, R_H_GeVinv))
report("Lambda_* = (rho_c R_H)^(1/3) [MeV]", 29.5, Lambda_star * 1e3, tol = 0.005)

cat("\nDelta N = 3 ln(Lambda_CP/Lambda_*), for Lambda_CP > Lambda_*:\n\n")

deltaN <- function(Lambda_CP_GeV) {
  ratio <- (Lambda_CP_GeV / Lambda_star)^3
  ifelse(ratio > 1, log(ratio), 0)
}

report("Delta N at Lambda_CP = 1 GeV",   10.6, deltaN(1),   tol = 0.005)
report("Delta N at Lambda_CP = 1e3 GeV", 31.3, deltaN(1e3), tol = 0.005)
report("Delta N at Lambda_CP = 1e8 GeV (paper rounds to 66)",
       65.8, deltaN(1e8), tol = 0.005)

cat("\nFull table (cp_wall_efolds_output.txt):\n")
cat(sprintf("  %12s %14s %14s\n", "Lambda_CP", "d/R_H needed", "extra e-folds"))
for (L in c(1e-3, 1e-2, 1, 1e3, 1e8, 1e12, 1e16)) {
  ratio <- (L / Lambda_star)^3
  dN <- deltaN(L)
  cat(sprintf("  %12.0e %14.2e %14.1f\n", L, ratio, dN))
}

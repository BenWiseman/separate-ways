#!/usr/bin/env Rscript
# neutrino_masses.R -- reproduces Sum(m_nu) and the m_betabeta range for the
# "one massless light neutrino" (BFT Z2) scenario, both orderings.
#
# SOURCE: beyond/push/P2_NUR_DM/fold_structure.py, section "4. light-neutrino
# observables with one exactly massless state". The oscillation inputs below are
# copied verbatim from that file's "ASSUMED oscillation inputs (NuFIT-6.0-era
# central values)" block, lines ~104-108, and cited here rather than re-derived:
#
#   dm21    = 7.49e-5  eV^2         (Delta m^2_21)
#   dm3l_NO = 2.513e-3 eV^2         (|Delta m^2_31|, normal ordering)
#   dm3l_IO = 2.484e-3 eV^2         (|Delta m^2_32|, inverted ordering)
#   sin^2(theta12) = 0.308
#   sin^2(theta13) = 0.02215 (NO), 0.02231 (IO)
#
# Formulas (fold_structure.py, same section), also displayed in
# pub/paper2/PAPER2_v3.md secs 5.2 and 5.6:
#
#   Normal ordering,   m1 = 0:  m2 = sqrt(dm21), m3 = sqrt(dm3l_NO)
#   Inverted ordering, m3 = 0:  NuFIT quotes Delta m^2_32, so
#                                m2 = sqrt(|dm3l_IO|),
#                                m1 = sqrt(|dm3l_IO| - dm21)
#   Sum m_nu = m1 + m2 + m3
#   m_betabeta = |sum_i U_ei^2 m_i|, extremal over Majorana phases so the two
#                bracketing values are | |c13^2 s12^2 m2 -/+ s13^2 m3| |  (NO)
#                and                   | (1-s12^2) c13^2 m1 -/+ s12^2 c13^2 m2 | (IO)
#                i.e. m_betabeta in [ |t2 - t3|, t2 + t3 ] with t2, t3 as below.
#
# Expected headline numbers (PAPER2_v3.md secs 5.2 and 5.6; fold_structure.py output):
#   normal,   m1=0: Sum m_nu = 58.8 meV,  m_betabeta in [1.5, 3.7] meV
#   inverted, m3=0: Sum m_nu = 98.9 meV, m_betabeta in [18.2, 48.2] meV
#
# No new physics: this is a direct port of fold_structure.py's arithmetic into R.

.args <- commandArgs(trailingOnly = FALSE)
.f <- sub("--file=", "", .args[grep("--file=", .args)])
.dir <- if (length(.f)) dirname(normalizePath(.f)) else "."
source(file.path(.dir, "helpers.R"))

section("neutrino_masses.R -- fold_structure.py sec 4, one massless light neutrino")

# ---------------------------------------------------------------------------
# Inputs, copied verbatim from beyond/push/P2_NUR_DM/fold_structure.py lines 104-108
# ---------------------------------------------------------------------------
dm21    <- 7.49e-5     # eV^2
dm3l_NO <- 2.513e-3    # eV^2
dm3l_IO <- 2.484e-3    # eV^2
s12sq   <- 0.308
s13sq_NO <- 0.02215
s13sq_IO <- 0.02231

cat(sprintf("\nInputs (fold_structure.py, NuFIT-6.0-era central values):\n"))
cat(sprintf("  dm21 = %.3e eV^2, |dm3l|(NO) = %.3e, |dm3l|(IO) = %.3e\n", dm21, dm3l_NO, dm3l_IO))
cat(sprintf("  sin^2(th12) = %s, sin^2(th13) = %s (NO), %s (IO)\n\n",
            s12sq, s13sq_NO, s13sq_IO))

# --- Normal ordering, m1 = 0 -------------------------------------------------
m1 <- 0; m2 <- sqrt(dm21); m3 <- sqrt(dm3l_NO)
sum_NO <- m1 + m2 + m3
c13sq_NO <- 1 - s13sq_NO
t2 <- s12sq * c13sq_NO * m2
t3 <- s13sq_NO * m3
mbb_NO <- c(abs(t2 - t3), t2 + t3)

cat(sprintf("NORMAL ordering, m_lightest = 0:\n"))
cat(sprintf("  m = (0, %.2f, %.2f) meV\n", m2 * 1e3, m3 * 1e3))

## Tolerances below are ABSOLUTE and matched to the number of decimal digits the
## source/paper actually prints (half a unit in the last printed digit, doubled for
## safety margin): all three headline values are printed to one decimal here, so
## the absolute rounding tolerance is +-0.05 meV. A relative-percent tolerance would be the wrong tool
## here -- it would flag correct bit-for-bit reproductions of a rounded display
## value as failures.
report("Sum m_nu, normal ordering [meV]", 58.8, sum_NO * 1e3, tol = 0.05, mode = "abs")
report("m_betabeta lower bound, NO [meV]", 1.5, mbb_NO[1] * 1e3, tol = 0.05, mode = "abs")
report("m_betabeta upper bound, NO [meV]", 3.7, mbb_NO[2] * 1e3, tol = 0.05, mode = "abs")

# --- Inverted ordering, m3 = 0 ------------------------------------------------
M3 <- 0; M1_ <- sqrt(dm3l_IO - dm21); M2_ <- sqrt(dm3l_IO)
sum_IO <- M1_ + M2_ + M3
c13sq_IO <- 1 - s13sq_IO
u1 <- (1 - s12sq) * c13sq_IO * M1_
u2 <- s12sq * c13sq_IO * M2_
mbb_IO <- c(abs(u1 - u2), u1 + u2)

cat(sprintf("\nINVERTED ordering, m_lightest = 0:\n"))
cat(sprintf("  m = (%.2f, %.2f, 0) meV\n", M1_ * 1e3, M2_ * 1e3))

report("Sum m_nu, inverted ordering [meV]", 98.9, sum_IO * 1e3, tol = 0.05, mode = "abs")
report("m_betabeta lower bound, IO [meV]", 18.2, mbb_IO[1] * 1e3, tol = 0.05, mode = "abs")
report("m_betabeta upper bound, IO [meV]", 48.2, mbb_IO[2] * 1e3, tol = 0.05, mode = "abs")

cat("\nCross-check against DESI DR2 BAO+CMB bound Sum m_nu < 64.2 meV (95%, [35] in PAPER2_v3.md):\n")
cat(sprintf("  NO: %.2f meV < 64.2 meV -> %s\n", sum_NO * 1e3, if (sum_NO*1e3 < 64.2) "ALIVE" else "EXCLUDED"))
cat(sprintf("  IO: %.2f meV < 64.2 meV -> %s\n", sum_IO * 1e3, if (sum_IO*1e3 < 64.2) "ALIVE" else "EXCLUDED"))

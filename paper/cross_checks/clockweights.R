#!/usr/bin/env Rscript
# clockweights.R -- reproduces the "how much time a field experiences" clock weights
# of TANGENTS_20260908.md section 22.
#
# SOURCE: TANGENTS_20260908.md sec 22 ("Ben's 'the graviton lives in the present'
# turned into the record's own numbers"), which builds on the exact per-mode weight
# w_n(A) = A^4(A^2-1)/(n^2(n^2-1)^2) of FIRST_TICK.md sec 2.4 and states the
# SUB-HORIZON asymptote used here directly:
#
#   w ~= (aH/k)^6 = (lambda / (2 pi lambda_H))^6 ,   lambda_H := c/H0
#
# ("checked: the exact/asymptote ratio is 1.0000 at k/aH = 3, 10, 100" -- every
# detector band below is enormously deeper into the sub-horizon regime than that,
# so the asymptote is exact for this purpose; the exact mode sum itself is not
# repeated here, only the asymptotic formula sec 22 gives).
#
# For a massive field, sec 22 quotes a (m/H)^2 law (from GRAVITY_CLOCK, "extrapolated
# to today's Hubble rate, not a new derivation") and gives hbar*H0 = 1.44e-33 eV.
#
# Quoted numbers (TANGENTS_20260908.md sec 22), reproduced below:
#   LIGO   (100 Hz)  : w = 1.8e-123
#   LISA   (1 mHz)   : w = 1.8e-93
#   PTA    (3 nHz)   : w = 2e-60
#   electron (m/H0)^2: 1.3e77
#   proton   (m/H0)^2: 4e83
#
# PRECISION NOTE, stated up front: the source quotes these to only 1-2 significant
# figures, and w depends on the input frequency (or mass) to a HIGH power (6th power
# for the GW-band numbers, 2nd power for the mass ones). A percent-level rounding of
# a "nominal" frequency like PTA's "3 nHz" (which is a representative order-of-magnitude
# PTA frequency, not a precision-quoted one) therefore produces a much larger swing in
# w. This script uses an order-of-magnitude ("oom", factor-based) comparison for these
# five numbers rather than the 0.5% default, and says so rather than picking a
# relative-percent tolerance that would look precise but isn't meaningful here.

.args <- commandArgs(trailingOnly = FALSE)
.f <- sub("--file=", "", .args[grep("--file=", .args)])
.dir <- if (length(.f)) dirname(normalizePath(.f)) else "."
source(file.path(.dir, "helpers.R"))

section("clockweights.R -- TANGENTS_20260908.md sec 22 graviton/fermion clock weights")

# ---------------------------------------------------------------------------
# Constants (all cited)
# ---------------------------------------------------------------------------
c_light <- 299792458            # m/s, exact SI definition
Mpc_m   <- 3.0857e22            # m per Mpc -- same conversion beyond/gut/calc_cpwall/cp_wall_efolds.py uses
H0_kms_Mpc <- 67.4              # km/s/Mpc, as stated in the record
H0 <- H0_kms_Mpc * 1000 / Mpc_m # s^-1

lambda_H <- c_light / H0
cat(sprintf("\nH0 = %.6e s^-1  (from %.1f km/s/Mpc)\n", H0, H0_kms_Mpc))
report("lambda_H = c/H0 [m]", 1.37e26, lambda_H, tol = 0.005,
       note = "TANGENTS sec22: 'lambda_H = c/H0 = 1.37e26 m'")

w_subhorizon <- function(freq_Hz) {
  lambda <- c_light / freq_Hz
  (lambda / (2 * pi * lambda_H))^6
}

cat("\nGraviton clock weight w = (aH/k)^6 = (lambda/2*pi*lambda_H)^6, evaluated today (a=1,H=H0):\n\n")

w_ligo <- w_subhorizon(100)          # LIGO band, 100 Hz
w_lisa <- w_subhorizon(1e-3)         # LISA band, 1 mHz
w_pta  <- w_subhorizon(3e-9)         # PTA band, 3 nHz

report("w(LIGO, 100 Hz)", 1.8e-123, w_ligo, tol = 1.5, mode = "oom")
report("w(LISA, 1 mHz)",  1.8e-93,  w_lisa, tol = 1.5, mode = "oom")
report("w(PTA, 3 nHz)",   2e-60,    w_pta,  tol = 1.5, mode = "oom",
       note = "largest quoted rounding of the five: 3 nHz is a nominal PTA frequency, and w~f^-6 amplifies it")

# ---------------------------------------------------------------------------
# Massive-fermion clock weight (m/H0)^2, hbar*H0 = 1.44e-33 eV (sec 22)
# ---------------------------------------------------------------------------
cat("\nMassive-fermion clock weight (m/H0)^2, via hbar*H0:\n\n")

hbar_eVs <- 6.582119569e-16      # eV s, CODATA (== dm_clock.py's GEV_S constant re-expressed in eV s)
hH0_eV <- hbar_eVs * H0
report("hbar*H0 [eV]", 1.44e-33, hH0_eV, tol = 0.01,
       note = "TANGENTS sec22: 'hbar H0 = 1.44e-33 eV'")

me_eV <- 0.51099895000e6         # electron mass, CODATA, eV
mp_eV <- 938.27208943e6          # proton mass, CODATA, eV

w_electron <- (me_eV / hH0_eV)^2
w_proton   <- (mp_eV / hH0_eV)^2

report("(m_e/H0)^2 [electron clock weight]", 1.3e77, w_electron, tol = 1.5, mode = "oom")
report("(m_p/H0)^2 [proton clock weight]",   4e83,   w_proton,   tol = 1.5, mode = "oom")

cat("\nSummary (all quoted to 1-2 sig figs in the source; 'oom' = order-of-magnitude,\n")
cat("factor-based check, since w depends on the 6th power of frequency or 2nd power of\n")
cat("mass, so 2-sig-fig input rounding is amplified far past 0.5%):\n")
cat(sprintf("  %-22s %14s %14s\n", "quantity", "expected", "reproduced"))
cat(sprintf("  %-22s %14s %14s\n", "w(LIGO)",  "1.8e-123", fmt_num(w_ligo)))
cat(sprintf("  %-22s %14s %14s\n", "w(LISA)",  "1.8e-93",  fmt_num(w_lisa)))
cat(sprintf("  %-22s %14s %14s\n", "w(PTA)",   "2e-60",    fmt_num(w_pta)))
cat(sprintf("  %-22s %14s %14s\n", "(m_e/H0)^2", "1.3e77", fmt_num(w_electron)))
cat(sprintf("  %-22s %14s %14s\n", "(m_p/H0)^2", "4e83",   fmt_num(w_proton)))

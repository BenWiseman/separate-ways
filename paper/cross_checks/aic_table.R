# aic_table.R -- parsimony accounting (AIC) for the framework against LCDM and w0waCDM, from the
# record's own geometry-only likelihood profiles (BAO + Pantheon+ + omega_cb + theta_* priors;
# calc/tangents/numass_mirror/numass_mirror_output.txt, lines cited). Base R only.
# AIC = chi2_min + 2k ;  Delta AIC = Delta chi2 - 2 Delta k  (positive = worse than the reference).
#
# Inputs (verbatim from numass_mirror_output.txt):
#   line 79: lcdm  chi2_min = 1406.859 at Sum m_nu = 0.0 meV
#   line 82: lcdm  Delta chi2 at Sum m_nu = 58.8 meV vs own best = +3.012 (= 1.74 sigma)
#   line 83: w0wa  chi2_min = 1397.739 at Sum m_nu = 0.0 meV
#   line 86: w0wa  Delta chi2 at Sum m_nu = 58.8 meV vs own best = +1.420
#   line 111: at Sum m_nu = 58.8 meV FIXED: w0wa - LCDM = -10.71 (2 params)
# Parameter bookkeeping relative to LCDM with Sum m_nu free (k0 parameters):
#   framework structure  : Sum m_nu fixed (58.8 meV), r = 0, Lambda (w = -1)  -> k0 - 1 on this data
#                          (r and w are not varied in the LCDM reference either, so only Sum m_nu counts here)
#   w0waCDM              : Sum m_nu free, w0, wa free                           -> k0 + 2
#   framework + tilt mech: additionally n_s fixed at 0.957888 (Turok-Boyle) -> k0 - 2, with the
#                          CMB tilt tension 4.83 sigma (Paper 2 sec 5.2, P-ACT-LB) as Delta chi2 = 4.83^2
source("helpers.R")
chi2_lcdm <- 1406.859; dchi2_numass_lcdm <- 3.012
chi2_w0wa <- 1397.739; dchi2_numass_w0wa <- 1.420
d_w0wa_vs_lcdm_fixed <- -10.71
# --- geometry-only, Sum m_nu free in the references ---
# AIC_fold - AIC_lcdm = (chi2_lcdm + 3.012 + 2(k0-1)) - (chi2_lcdm + 2 k0) = 3.012 - 2
dAIC_fold_vs_lcdm  <- dchi2_numass_lcdm - 2
dAIC_w0wa_vs_lcdm  <- (chi2_w0wa - chi2_lcdm) + 4                       # Dk = +2
dAIC_fold_vs_w0wa  <- dAIC_fold_vs_lcdm - dAIC_w0wa_vs_lcdm
# cross-check with line 111 (both at 58.8 meV fixed): w0wa(fixed) - lcdm(fixed) = -10.71 with Dk = +2
dAIC_w0wafixed_vs_fold <- d_w0wa_vs_lcdm_fixed + 4
# --- adding the Turok-Boyle tilt mechanism (n_s fixed): CMB tilt tension only, sec 5.2 ---
tilt_sigma <- 4.83; dAIC_mech_vs_lcdm_tilt <- tilt_sigma^2 - 2
cat(sprintf("Delta AIC (framework structure vs LCDM, geometry-only)      = %+.2f   [Dchi2 %+.3f, Dk = -1]\n", dAIC_fold_vs_lcdm, dchi2_numass_lcdm))
cat(sprintf("Delta AIC (w0waCDM vs LCDM, geometry-only)                  = %+.2f   [Dchi2 %+.3f, Dk = +2]\n", dAIC_w0wa_vs_lcdm, chi2_w0wa-chi2_lcdm))
cat(sprintf("Delta AIC (framework structure vs w0waCDM, geometry-only)   = %+.2f\n", dAIC_fold_vs_w0wa))
cat(sprintf("cross-check at Sum m_nu fixed for both: w0wa vs framework   = %+.2f   [line 111: -10.71, Dk = +2]\n", dAIC_w0wafixed_vs_fold))
cat(sprintf("Delta AIC (adding the tilt mechanism, n_s fixed; CMB tilt)  = %+.1f   [%.2f sigma -> Dchi2 %.1f, Dk = -1 more]\n", dAIC_mech_vs_lcdm_tilt, tilt_sigma, tilt_sigma^2))
# --- BIC on the same data: n = 1580 Pantheon+ SNe (zHD > 0.01, non-calibrator; numass_mirror.py l.147)
#     + 13 DESI DR2 BAO points (boltzmann_twoply._desi_vectors) + 2 compressed priors = 1595 ---
n_data <- 1580 + 13 + 2; pen <- log(n_data)
dBIC_fold_vs_lcdm <- dchi2_numass_lcdm - pen
dBIC_w0wa_vs_lcdm <- (chi2_w0wa - chi2_lcdm) + 2*pen
dBIC_fold_vs_w0wa <- dBIC_fold_vs_lcdm - dBIC_w0wa_vs_lcdm
cat(sprintf("\nBIC, n = %d, ln n = %.3f per parameter:\n", n_data, pen))
cat(sprintf("Delta BIC (framework structure vs LCDM)   = %+.2f\n", dBIC_fold_vs_lcdm))
cat(sprintf("Delta BIC (w0waCDM vs LCDM)               = %+.2f\n", dBIC_w0wa_vs_lcdm))
cat(sprintf("Delta BIC (framework structure vs w0waCDM) = %+.2f\n", dBIC_fold_vs_w0wa))
# --- Published, CMB-inclusive block (DESI DR2; verbatim in working note THEORY_COMPARE_20260912.md) ---
#   Elbers et al. 2503.14744 Table 4, DESI+CMB: fixing Sum m_nu at the oscillation floor costs
#   Delta chi2_MAP = 7.2 (goodness-of-fit-loss metric, 2.7 sigma).
#   DESI DR2 2503.14738 Table 6, w0waCDM vs LCDM Delta chi2_MAP: DESI+CMB -12.5 (3.1s), +Pantheon+ -10.7 (2.8s),
#   +Union3 -17.4 (3.8s), +DESY5 -21.0 (4.2s).
# The effective number of independent data points in a CMB-inclusive likelihood is not a simple count, so
# the BIC verdict is stated as the BREAK-EVEN n: BIC prefers the extra parameters only if n_eff < exp(-Dchi2/(2 Dk)).
dchi2_pub_fold <- 7.2
w0wa_pub <- c("DESI+CMB"=-12.5, "DESI+CMB+Pantheon+"=-10.7, "DESI+CMB+Union3"=-17.4, "DESI+CMB+DESY5"=-21.0)
cat("\nPublished CMB-inclusive block:\n")
cat(sprintf("framework vs LCDM: Dchi2 = +%.1f, Dk = -1 -> Delta AIC = %+.1f ; BIC prefers the framework iff n_eff > %.0f\n",
            dchi2_pub_fold, dchi2_pub_fold - 2, exp(dchi2_pub_fold)))
for (nm in names(w0wa_pub)) {
  d <- w0wa_pub[[nm]]
  cat(sprintf("w0waCDM vs LCDM, %-22s Dchi2 = %+.1f, Dk = +2 -> Delta AIC = %+.1f ; BIC prefers w0wa iff n_eff < %.0f\n",
              nm, d, d + 4, exp(-d/2)))
}
cat("\nchecks:\n")
report("Delta AIC framework vs LCDM (paper: 'a tie, +1')", expected = 1.0, reproduced = dAIC_fold_vs_lcdm, tol = 0.02, mode = "rel",
       note = "3.012 - 2 = 1.012; paper quotes +1")
report("Delta AIC w0waCDM vs LCDM, geometry-only (paper: '-5')", expected = -5.1, reproduced = dAIC_w0wa_vs_lcdm, tol = 0.01, mode = "rel",
       note = "-9.12 + 4 = -5.12; DESI's own CMB-inclusive values are stronger (see THEORY_COMPARE note)")
report("Delta AIC tilt mechanism vs LCDM (paper: '+21')", expected = 21.3, reproduced = dAIC_mech_vs_lcdm_tilt, tol = 0.01, mode = "rel",
       note = "4.83^2 - 2 = 21.3")
report("Delta BIC framework vs LCDM (paper: '-4.4')", expected = -4.4, reproduced = dBIC_fold_vs_lcdm, tol = 0.02, mode = "rel",
       note = "3.012 - ln(1595)")
report("Delta BIC w0waCDM vs LCDM (paper: '+5.6')", expected = 5.6, reproduced = dBIC_w0wa_vs_lcdm, tol = 0.02, mode = "rel",
       note = "-9.12 + 2 ln(1595)")
report("break-even n, framework vs LCDM, published (paper: '1300')", expected = 1340, reproduced = exp(dchi2_pub_fold), tol = 0.01, mode = "rel",
       note = "exp(7.2) = 1339.4")
report("break-even n, w0wa vs LCDM, DESI+CMB+Pantheon+ (paper: '210')", expected = 211, reproduced = exp(10.7/2), tol = 0.01, mode = "rel",
       note = "exp(5.35) = 210.6")
report("break-even n, w0wa vs LCDM, DESI+CMB+DESY5 (paper: '36000')", expected = 36316, reproduced = exp(21.0/2), tol = 0.01, mode = "rel",
       note = "exp(10.5) = 36315.5")

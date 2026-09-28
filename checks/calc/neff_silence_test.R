#!/usr/bin/env Rscript
# neff_silence_test.R
#
# TWO THINGS. First, my own partner_sheet_test.R is wrong in its statistics and
# this corrects it. Second, an adversarial pass pointed at a far better test that
# uses EXISTING data, and this prices it.
#
# THE ERROR IN partner_sheet_test.R: I computed the events needed as
# sigma_mean = resolution / sqrt(N), which is right for random per-event
# resolution and WRONG when the limiting uncertainty is a systematic shift of the
# absolute energy scale common to every event. A common-mode systematic does not
# average down. At EeV the absolute cascade energy scale carries contributions
# from the inelasticity distribution at x ~ 1e-6, the ice absorption model, PMT
# saturation and the background shape; the adversarial estimate is 15-25 per cent
# in total, against a signal separation of 7.18 per cent.
#
# THE BETTER TEST: if the partner sheet is dynamically silent, it contributes no
# radiation density at ANY epoch, including recombination, where N_eff is
# measured directly.

cat("=== 1. correcting my own error: a systematic does not average down ===\n")
frac_sep  <- 2^(1/10) - 1        # 7.18 per cent, the signal
res_rand  <- 0.30                # per-event random resolution
syst_lo   <- 0.15; syst_hi <- 0.25   # absolute energy-scale systematic

cat(sprintf("  signal separation                 : %.2f per cent\n", 100*frac_sep))
cat(sprintf("  per-event random resolution       : %.0f per cent  (averages as 1/sqrt(N))\n", 100*res_rand))
cat(sprintf("  absolute energy-scale systematic  : %.0f-%.0f per cent  (does NOT average)\n",
            100*syst_lo, 100*syst_hi))
for (N in c(100, 160, 1000, 1e5)) {
  tot_lo <- sqrt((res_rand/sqrt(N))^2 + syst_lo^2)
  tot_hi <- sqrt((res_rand/sqrt(N))^2 + syst_hi^2)
  cat(sprintf("  N = %-7s combined uncertainty %.1f-%.1f per cent -> significance %.2f-%.2f sigma\n",
              format(N, big.mark=","), 100*tot_lo, 100*tot_hi,
              frac_sep/tot_hi, frac_sep/tot_lo))
}
cat("\n  So the neutrino-line test is systematics-limited, not statistics-limited.\n")
cat(sprintf("  It needs absolute calibration better than %.1f per cent to reach 1 sigma at all,\n",
            100*frac_sep))
cat(sprintf("  and better than %.2f per cent for 3 sigma. No existing detector is close.\n",
            100*frac_sep/3))

cat("\n=== 2. the test that uses data we already have ===\n")
# radiation density relative to photons: rho_rad = rho_gamma (1 + 0.2271 N_eff)
coef <- 0.2271
Neff_sm <- 3.044
rho_over_gamma <- function(Neff) 1 + coef*Neff
cat(sprintf("  Standard Model radiation density  : %.4f rho_gamma at N_eff = %.3f\n",
            rho_over_gamma(Neff_sm), Neff_sm))
# a loud partner sheet contributes an identical radiation budget that we do not
# see as photons, so it appears entirely as extra N_eff
extra <- rho_over_gamma(Neff_sm)          # in units of rho_gamma
dNeff <- extra/coef
Neff_loud <- Neff_sm + dNeff
cat(sprintf("  a loud partner sheet adds         : %.4f rho_gamma = %.2f in N_eff units\n",
            extra, dNeff))
cat(sprintf("  predicted N_eff if the sheet is loud: %.2f\n", Neff_loud))

Neff_obs <- 2.99; Neff_err <- 0.17       # Planck
cat(sprintf("\n  measured (Planck)                 : %.2f +/- %.2f\n", Neff_obs, Neff_err))
cat(sprintf("  silent prediction                 : %.3f  -> %.2f sigma from the measurement\n",
            Neff_sm, abs(Neff_sm-Neff_obs)/Neff_err))
cat(sprintf("  loud prediction                   : %.2f  -> %.1f sigma from the measurement\n",
            Neff_loud, abs(Neff_loud-Neff_obs)/Neff_err))

cat("\n=== 3. validation: the check must be able to fail ===\n")
stopifnot(abs(Neff_sm - Neff_obs)/Neff_err < 1)      # silent must be consistent
stopifnot(abs(Neff_loud - Neff_obs)/Neff_err > 10)   # loud must be excluded
bad <- tryCatch({ stopifnot(abs(Neff_sm-Neff_obs)/Neff_err > 10); TRUE }, error=function(e) FALSE)
cat(sprintf("  silent consistent, loud excluded, and the reversed assertion fails: %s\n",
            if (!bad) "yes" else "NO - BLIND"))

cat(sprintf("
=== 4. flatly ===

  The neutrino line is the wrong place to test the partner sheet's silence and my
  earlier script got there by treating a systematic as if it averaged down. The
  separation is %.1f per cent of the line energy and the absolute energy scale at
  these energies is uncertain at %.0f-%.0f per cent, so the significance saturates
  below one sigma no matter how many events arrive.

  The right test already exists. If the partner sheet is dynamically silent it
  contributes no radiation at any epoch, and at recombination that quantity is
  measured directly as N_eff. A loud sheet carrying the same content adds its
  whole radiation budget without adding photons we can see, which is exactly what
  N_eff counts: it predicts %.1f against a measured %.2f +/- %.2f, excluded at
  %.0f sigma. The silent reading predicts %.3f and sits %.1f sigma from the
  measurement.

  So the commitment this paper makes about what the second sheet IS has already
  been tested, and it passed. That is a stronger thing to be able to say than a
  future measurement, and it belongs at the front of the paper rather than in a
  caveat.

  THE CAVEAT THAT MUST TRAVEL WITH IT. This tests silence AT RECOMBINATION. The
  cosmology paper's g_* commitment is about the bang. If silence is a property at
  all times the test applies directly; if the sheet could be loud at the bang and
  silent later, the framework owes an account of when silence turns on, and this
  test constrains only the late half. Say which, rather than letting a referee ask.
", 100*frac_sep, 100*syst_lo, 100*syst_hi,
   Neff_loud, Neff_obs, Neff_err, abs(Neff_loud-Neff_obs)/Neff_err,
   Neff_sm, abs(Neff_sm-Neff_obs)/Neff_err))

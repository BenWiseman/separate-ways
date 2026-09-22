# t_dec ~ M_1^(2/3) is the strongest single relation in the paper, and the DESI corner sits
# beside it. The paper says the LCDM bound on Sigma m_nu is "nine per cent above" the
# 58.8 meV floor and that the pair is testable at the next release. Nine per cent of what,
# and how much tightening kills it? The paper states the margin but never converts it into
# a survival condition, which is the form a referee can act on.

floor <- 58.8                     # meV, exact-stabilisation floor with normal ordering
bound <- 64.2   # meV. CORRECTED 2026-09-21: was 64.0, the rounded form, which put
                # 5.2 meV and 8.1 per cent into the paper against 5.4 and 8.4
                # everywhere else. A later check caught the drift. VERIFIED against the manuscript: Elbers et al., DESI DR2 BAO + CMB,
                # quoted in 3.3 as "the same data give Sigma m_nu < 64.2". Not inferred from the
                # "nine per cent" phrasing.
cat(sprintf("=== 1. the margin as it stands\n\n   floor %.1f meV, bound %.1f meV, gap %.1f meV = %.1f%% of the floor\n",
    floor, bound, bound-floor, 100*(bound-floor)/floor))

cat("\n=== 2. the survival condition, which is what the paper does not state\n\n")
cat("   The model dies when the bound falls BELOW the floor. So it needs the bound to\n")
cat("   tighten by more than the gap. Expressed as the fractional improvement required:\n\n")
cat("        bound falls to (meV)   fractional tightening    model status\n")
for (b in c(64.0, 62.0, 60.0, 58.8, 57.0, 54.0)) 
  cat(sprintf("   %20.1f %22.1f%% %16s\n", b, 100*(bound-b)/bound,
      if (b > floor) "survives" else if (abs(b-floor) < 1e-9) "exactly at the floor" else "REFUTED"))
cat(sprintf("\n   So the whole question is an %.1f%% tightening of one published number.\n",
    100*(bound-floor)/bound))

cat("\n=== 3. is that a plausible amount for one data release?\n\n")
cat("   Neutrino-mass bounds scale roughly as the inverse square root of survey volume\n")
cat("   for a fixed systematic floor, so a bound improves by about (V_new/V_old)^-1/2.\n\n")
cat("        volume ratio    expected bound (meV)    crosses the floor?\n")
for (r in c(1.00, 1.05, 1.10, 1.19, 1.30, 1.50)) {
  b <- bound/sqrt(r)
  cat(sprintf("   %14.1f %22.1f %20s\n", r, b, if (b < floor) "YES" else "no"))
}
cat(sprintf("\n   The floor is crossed once the volume grows by a factor of %.2f.\n", (bound/floor)^2))
cat("   That is a modest multiple, not a generational one, which is why the paper is\n")
cat("   right that this is testable at the next release rather than eventually.\n")

cat("\n=== 4. what this adds and what it does not\n\n")
cat("  ADDS a survival condition in the form a referee can use: the model requires the\n")
cat("  published bound NOT to tighten by more than 8.4 per cent, equivalently requires the\n")
cat("  effective survey volume not to grow by more than about 1.19 at fixed systematics.\n")
cat("  That is a sharper statement than 'within nine per cent' and it is falsifiable on a\n")
cat("  named timescale.\n")
cat("  DOES NOT add a forecast. The volume scaling is a rule of thumb, the real bound\n")
cat("  depends on the dataset combination, the priors and the systematic floor, and this\n")
cat("  computes none of those. It converts the margin the paper already quotes into a\n")
cat("  condition; it does not predict when the condition will be met.\n")

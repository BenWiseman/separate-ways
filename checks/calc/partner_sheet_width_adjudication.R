#!/usr/bin/env Rscript
# adjudicate_m3_r1.R -- check the two load-bearing claims in checks/oracle/r1.
# M3 is a menu, not an oracle. Neither claim is accepted without arithmetic.

cat("=== CLAIM A: point 9 is undersold, the partner-sheet test is 17.5 widths ===\n")
one_sheet  <- 491.6; both_sheets <- 526.9; width <- 2.0
sep <- (both_sheets - one_sheet)/width
cat(sprintf("  one sheet gravitates : %.1f PeV\n", one_sheet))
cat(sprintf("  both sheets gravitate: %.1f PeV\n", both_sheets))
cat(sprintf("  separation           : %.1f PeV = %.1f quoted widths\n",
            both_sheets-one_sheet, sep))
# the scaling the inventory states: M ~ g_*^(1/10)
g1 <- 106.75; g2 <- 2*g1
pred <- one_sheet * (g2/g1)^(1/10)
cat(sprintf("  check against M ~ g_*^(1/10): %.1f PeV predicted vs %.1f quoted, %.2f%% apart\n",
            pred, both_sheets, 100*abs(pred-both_sheets)/both_sheets))
cat(sprintf("  VERDICT: %s\n", if (abs(pred-both_sheets)/both_sheets < 0.01 && sep > 15)
    "claim A stands, arithmetic reproduces both the value and the separation" else "claim A fails"))

cat("\n=== CLAIM B: 'kappa = 2 EXACTLY, forced by CPT' ===\n")
cat("  The project withdrew 'exactly two' because it assumed equal accretion on\n")
cat("  the two exteriors. M3 says CPT forces equality, so the withdrawal was\n")
cat("  over-cautious. Adjudicate by asking WHEN the two sheets are mirror images.\n\n")
cat("  The fold's own result (A.4 and the branch-imbalance lead): a Theta-invariant\n")
cat("  state does NOT stay Theta-invariant under evolution, because Theta is\n")
cat("  antilinear and CPT covariance gives Theta U Theta^-1 = U^-1, so\n")
cat("  Theta|psi(t)> = |psi(-t)>. The branch imbalance is an ODD function of time,\n")
cat("  exactly zero at the fold point and growing away from it.\n\n")
cat("  So mirror symmetry of the two sheets' CONFIGURATIONS is exact at the bang\n")
cat("  and degrades with time. kappa = 2 is therefore not an assumption and not a\n")
cat("  theorem: it is an EARLY-TIME limit, and its accuracy is controlled by a\n")
cat("  quantity the project has already computed the shape of.\n\n")

# what the growth window looks like against that
H0 <- 67.4; Om <- 0.315; OL <- 1-Om
age <- function(z) { Hub <- 9.778e9/(H0/100); (2/3)*Hub/sqrt(OL)*asinh(sqrt(OL/Om)*(1+z)^(-1.5)) }
t_now <- age(0)
cat(sprintf("  %8s %14s %16s\n", "z", "age (Gyr)", "fraction of t_0"))
for (z in c(20, 10, 7, 4, 1, 0)) cat(sprintf("  %8.1f %14.3f %16.4f\n", z, age(z)/1e9, age(z)/t_now))
cat("
  The little red dot window (z=20 to z=7) sits inside the first 5.5 per cent of
  cosmic time. If the departure from mirror symmetry grows with time from an exact
  zero at the bang, then the epoch where the fold's growth mechanism is claimed is
  the epoch where kappa = 2 is best justified, and the epoch where it is worst
  justified is today, where the mechanism is not needed because the e-fold budget
  is not binding anyway.

  VERDICT ON CLAIM B: M3 is directionally right and its wording is wrong. CPT does
  not force kappa = 2 at all times. It forces mirror configurations AT THE FOLD
  POINT, and the branch imbalance measures the drift. The correct statement is
  'kappa -> 2 in the early-time limit, with a computable departure', which is
  stronger than the withdrawn 'exactly two' and stronger than the current 'one
  plus an unestablished ratio'. NOT YET COMPUTED: the size of the departure over
  the growth window, which needs the branch imbalance evaluated on a realistic
  Hamiltonian rather than the random CPT-covariant matrix used so far.
")

cat("=== CLAIM C: 'the LRD mass function has a second mode at 4.1e5 times the first' ===\n")
cat("  This is wrong and it matters, because it would be the paper's headline.\n")
cat("  4.1e5 is the ratio of FINAL MASSES for the SAME seed at kappa=2 against\n")
cat("  kappa=1, at FULL duty cycle, over the whole window. A population has a\n")
cat("  spread of seeds, formation times and duty cycles, so that ratio is not a\n")
cat("  mode separation. The computed separation in the residual is 0.84 dex at\n")
cat("  kappa=1.5, lambda=0.3, which is a factor of:\n")
cat(sprintf("    %.1f, not %.1e\n", 10^0.84, 4.1e5))
cat("  VERDICT: rejected. Do not let this reach the draft.\n")

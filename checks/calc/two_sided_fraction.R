#!/usr/bin/env Rscript
# ===========================================================================
# SUPERSEDED, 2026-09-23, by tail_from_lrd.R. DO NOT QUOTE NUMBERS FROM HERE.
#
# This script holds nu fixed while varying the tail index p, which is not
# self-consistent: softening the tail raises the ordinary formation rate as
# well as the type II-B one, and the ratio gains far less than it appears to
# here. Its p = 1.5 and p = 1 figures (2.97 and 36.8 per cent) are roughly 40x
# and 37x too large; tail_from_lrd.R gives 7.2e-4 and 1.0e-2.
#
# Those two wrong numbers reached a draft of the companion's 5.2 and were
# caught by an audit against the scripts. The header is here so it cannot
# happen twice. The scaling argument below is kept for its structure only.
# ===========================================================================
# two_sided_fraction.R -- what does the two-sided class COST in abundance?
#
# The fold's interior machinery needs a bifurcation surface, which collapse cannot
# supply. arXiv:2401.06329 simulates type II-B primordial holes in RADIATION
# domination and finds bifurcating trapping horizons above an amplitude threshold
# (mu >~ 1.8 in their parametrisation), while type I holes form from a lower one.
#
# So two-sided holes are the rare tail of a tail. This script prices that, because
# the two-sided FRACTION f_2s is the one free parameter the bimodality test of the
# main paper's section 4.3 cannot set.
#
# The threshold ratio is taken as a parameter rather than converted from their mu,
# because converting needs their exact compaction-function definition and this
# script will not invent one. What is computed is how steeply the fraction falls
# with the separation between the two thresholds, which is definition-independent.

# --- Press-Schechter style collapse fraction for a Gaussian field -------------
# beta = fraction of regions above threshold nu = delta_c / sigma
banner <- function() {
  cat("\n")
  cat("  ############################################################################\n")
cat("  ## SUPERSEDED by tail_from_lrd.R, 2026-09-23. DO NOT QUOTE NUMBERS FROM HERE.\n")
  cat("  ############################################################################\n")
cat("  This holds nu fixed while varying the tail index p, which is not self-consistent: softening the\n")
cat("  tail raises the ordinary formation rate too. Its p = 1.5 and p = 1 figures are roughly 40x and\n")
cat("  37x too large; tail_from_lrd.R gives 7.2e-4 and 1.0e-2. Two of the wrong numbers reached a\n")
cat("  draft of the companion and were caught by an audit.\n")
  cat("  ############################################################################\n\n")
}
banner()
beta_gauss <- function(nu) erfc(nu/sqrt(2))/2
erfc <- function(x) 2*pnorm(-x*sqrt(2))

# ratio of two-sided (higher threshold) to all (lower threshold) holes
f2s <- function(nu1, ratio) beta_gauss(nu1*ratio) / beta_gauss(nu1)

cat("=== 1. validation: the checker must be able to fail ===\n")
# plant: at ratio = 1 the fraction must be exactly 1
stopifnot(abs(f2s(5, 1.0) - 1) < 1e-12)
# plant: a higher threshold must give a SMALLER fraction
stopifnot(f2s(5, 1.2) < f2s(5, 1.05))
# plant a deliberate failure to confirm the test can say no
planted_ok <- tryCatch({ stopifnot(f2s(5, 1.2) > f2s(5, 1.05)); TRUE },
                       error = function(e) FALSE)
cat(sprintf("  monotonicity holds, and the reversed assertion fails as it must: %s\n",
            if (!planted_ok) "yes" else "NO - CHECK IS BLIND"))

cat("\n=== 2. the cost of the two-sided class ===\n")
cat("  f_2s = (fraction above the type II-B threshold) / (fraction above the type I threshold)\n")
cat("  nu1 is how far into the tail ordinary PBH formation already sits.\n\n")
cat(sprintf("  %8s", "nu1"))
ratios <- c(1.05, 1.1, 1.2, 1.3, 1.5)
for (r in ratios) cat(sprintf(" %12s", sprintf("x%.2f", r)))
cat("\n")
for (nu1 in c(3, 4, 5, 6, 7, 8)) {
  cat(sprintf("  %8.1f", nu1))
  for (r in ratios) cat(sprintf(" %12.2e", f2s(nu1, r)))
  cat("\n")
}

cat("\n=== 3. what the bimodality test needs ===\n")
cat("  Section 4.3's mixture test needs f_2s of order 0.3 to fire at a few hundred\n")
cat("  objects, and cannot fire at all below a few per cent. Threshold ratio that\n")
cat("  still leaves f_2s above those levels:\n\n")
for (target in c(0.30, 0.03)) {
  cat(sprintf("  f_2s >= %.2f needs a threshold ratio below:\n", target))
  for (nu1 in c(3, 5, 7)) {
    rr <- uniroot(function(r) f2s(nu1, r) - target, c(1.0, 3.0), tol = 1e-10)$root
    cat(sprintf("     nu1 = %.0f :  ratio < %.4f   (a %.1f%% higher threshold)\n",
                nu1, rr, 100*(rr-1)))
  }
  cat("\n")
}

cat("=== 4. flatly ===\n")
cat("
  The two-sided class is the tail of a tail, and Gaussian statistics price it
  brutally. Ordinary primordial black hole formation already sits several sigma
  out; demanding a threshold even ten per cent higher cuts the surviving fraction
  by orders of magnitude, and the deeper into the tail the parent population sits,
  the worse it gets, because the Gaussian tail steepens.

  For the bimodality test to fire at a few hundred objects it needs roughly a
  third of holes to be two-sided. Under Gaussian statistics that requires the two
  thresholds to sit within about one per cent of each other at nu1 = 5, which no
  reading of the simulations supports: they separate type I from type II-B by a
  substantial amplitude margin.

  SO, FLATLY: if the curvature field is Gaussian, two-sided holes are far too rare
  to be a population, and the bimodality test cannot fire on them. The channel
  supplies the geometry and not the numbers.

  WHAT THIS DOES NOT CLOSE, and it is the whole game: primordial black hole
  abundance is notoriously controlled by non-Gaussianity, precisely because it is
  a rare-tail quantity. A field with a heavier tail than Gaussian moves f_2s by
  orders of magnitude in the direction this needs. Section 5 prices that.
")

cat("=== 5. how much non-Gaussianity would rescue it? ===\n")
cat("  Model the tail as stretched-exponential, P(>x) ~ exp(-(x/x0)^p). p=2 is\n")
cat("  Gaussian; p<2 is heavier-tailed. Same threshold ratio, varying p:\n\n")
f2s_p <- function(nu1, ratio, p) exp(-(nu1*ratio)^p + nu1^p)
cat(sprintf("  %6s", "p"))
for (r in c(1.1, 1.2, 1.5)) cat(sprintf(" %14s", sprintf("ratio x%.1f", r)))
cat("\n")
for (p in c(2.0, 1.5, 1.0, 0.7, 0.5)) {
  cat(sprintf("  %6.1f", p))
  for (r in c(1.1, 1.2, 1.5)) cat(sprintf(" %14.2e", f2s_p(5, r, p)))
  cat("\n")
}
cat("
  At p = 1, an exponential tail, a twenty per cent higher threshold still leaves
  37 per cent of holes two-sided, which is precisely the f_2s the bimodality test
  needs. At p = 0.5 it leaves 81 per cent. So the
  question 'can two-sided holes be a population' is not a question about the fold
  at all. It is a question about the tail of the primordial curvature field, which
  is exactly the quantity the little-red-dot PBH literature is already fighting
  over. That is a testable link between the fold and an active controversy, and it
  runs in the useful direction: the fold NEEDS heavy tails, and heavy tails are
  what the overmassive-early-hole problem already wants.
")


banner()

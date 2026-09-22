# Does the bootstrap likelihood-ratio test rescue the two-sided black-hole test?
#
# deconvolved_mixture.R established two things. The NAIVE mixture is invalid once per-object
# errors vary: it reports two classes on one-class data in 89 per cent of draws at n=300. The
# DECONVOLVED fit, with the errors carried in the likelihood, is valid but weak: null 0.000 and
# power only 0.38 by n=3000. That left 4.3 saying the test is real and dearer than the idealised
# figure.
#
# BIC is not the only way to choose between one component and two, and in latent-class work it is
# not the best one. The bootstrap likelihood-ratio test -- fit both models, then get the null
# distribution of the LRT statistic by parametric bootstrap from the fitted ONE-component model
# rather than assuming chi-squared, which does not hold when a mixture weight sits on the boundary
# -- is the standard alternative in psychometrics (McLachlan 1987; Nylund, Asparouhov and Muthen
# 2007 report it beating BIC for choosing the number of classes). This file asks whether that
# carries across to a mixture with KNOWN heteroscedastic measurement errors, which is our case.
#
# This is the borrowable calculation the outreach to a psychometrician turns on, so it is run
# rather than cited: if BLRT does not beat BIC here, the ask is weaker and should be stated that
# way.
#
# Base R only, per the repo rule: no package is loaded. Parallelism is at the shell level
# (chunks run as separate Rscript processes), not via a library.
#
# It also settles a provenance question in 4.3. The claim that the test "fails at every size
# tried when the far side supplies a fifth" came from two_sided_bimodality_power.R, which is the
# superseded NAIVE analysis, and sits in a passage that now reads as if it were about the
# deconvolved estimator. Rather than argue which setting is harder, RATIO below runs the fifth
# case through the valid estimator directly.
#
# Usage:  Rscript blrt_vs_bic.R [n] [mode power|null] [seed] [ndraws] [B] [ratio]
#         ratio is the far side's accretion rate as a fraction of ours: 0.5 = half (kappa=1.5),
#         0.2 = a fifth (kappa=1.2). With no arguments it runs a small demonstration.

args <- commandArgs(trailingOnly = TRUE)
N      <- if (length(args) >= 1) as.integer(args[1]) else 300
MODE   <- if (length(args) >= 2) args[2] else "power"
SEED   <- if (length(args) >= 3) as.integer(args[3]) else 1
NDRAW  <- if (length(args) >= 4) as.integer(args[4]) else 5
B      <- if (length(args) >= 5) as.integer(args[5]) else 99
RATIO  <- if (length(args) >= 6) as.numeric(args[6]) else 0.5
# CONC = 0 keeps every class-2 object at the mean excess, which is what deconvolved_mixture.R
# and the first version of this file did. CONC > 0 draws a PER-OBJECT duty cycle from
# Beta(mean 0.3, concentration CONC), which is what two_sided_bimodality_power.R does and what
# the physics is: the excess is exp((kappa-1) lambda N1), so a low-duty-cycle object sits at
# zero excess on top of class 1. At CONC = 4, 27 per cent of class 2 sits below one intrinsic
# scatter. Ignoring that spread overstates the power, which is why this switch exists.
CONC   <- if (length(args) >= 7) as.numeric(args[7]) else 0
set.seed(20260921 + SEED)

# far side at RATIO of our rate, <lambda>=0.3, 12.93 e-folds: 0.5 -> 0.842 dex, 0.2 -> 0.337 dex
SEP  <- RATIO * 0.3 * 12.93 / log(10)
SINT <- 0.40                          # intrinsic scatter, dex
F2S  <- 0.30                          # two-sided fraction
het  <- function(n) 0.20 + rexp(n, rate = 1/0.40)   # per-object errors, a long tail of poor ones

# --- likelihoods with the per-object errors KNOWN (the deconvolved estimator)
nll2 <- function(par, x, sig) {
  p <- 1/(1+exp(-par[1])); m1 <- par[2]; m2 <- par[3]
  v1 <- exp(par[4])^2 + sig^2; v2 <- exp(par[5])^2 + sig^2
  -sum(log(pmax(p*dnorm(x, m1, sqrt(v1)) + (1-p)*dnorm(x, m2, sqrt(v2)), 1e-300)))
}
nll1 <- function(par, x, sig)
  -sum(dnorm(x, par[1], sqrt(exp(par[2])^2 + sig^2), log = TRUE))

fit1 <- function(x, sig) {
  o <- optim(c(mean(x), log(max(sd(x)*0.9, 1e-3))), nll1, x = x, sig = sig,
             method = "Nelder-Mead")
  list(val = o$value, par = o$par)
}
fit2 <- function(x, sig, tries = 4) {
  sd0 <- sd(x); best <- NULL
  for (t in seq_len(tries)) {
    q  <- as.numeric(quantile(x, c(0.3, 0.8))) + rnorm(2, 0, 0.1*sd0)
    st <- c(0, q[1], q[2], log(sd0*0.6), log(sd0*0.6))
    o  <- try(optim(st, nll2, x = x, sig = sig, method = "Nelder-Mead",
                    control = list(maxit = 2000)), silent = TRUE)
    if (!inherits(o, "try-error") && (is.null(best) || o$value < best$value)) best <- o
  }
  list(val = best$value, par = best$par)
}

# --- one dataset
N1 <- 12.93                                   # one-sided e-folds, z=20 to z=7
excess <- function(n2) {
  if (CONC <= 0) return(rep(SEP, n2))         # every object at the mean duty cycle
  lam <- rbeta(n2, 0.3*CONC, 0.7*CONC)        # per-object duty cycle, mean 0.3
  RATIO * lam * N1 / log(10)                  # dex; low lambda sits at zero excess
}
draw <- function(n, two_class) {
  sig <- het(n)
  if (two_class) { n2 <- rbinom(1, n, F2S); mu <- c(rep(0, n - n2), excess(n2)) }
  else           { mu <- rep(0, n) }
  list(x = rnorm(n, mu, sqrt(SINT^2 + sig^2)), sig = sig)
}

# --- the two decision rules on one dataset
decide <- function(x, sig, B) {
  f1 <- fit1(x, sig); f2 <- fit2(x, sig)
  n  <- length(x)
  bic1 <- 2*f1$val + 2*log(n); bic2 <- 2*f2$val + 5*log(n)
  lrt  <- 2*(f1$val - f2$val)
  # parametric bootstrap under the FITTED one-component model, reusing the same known errors
  mu0 <- f1$par[1]; s0 <- exp(f1$par[2])
  nulls <- numeric(B)
  for (b in seq_len(B)) {
    xb <- rnorm(n, mu0, sqrt(s0^2 + sig^2))
    g1 <- fit1(xb, sig); g2 <- fit2(xb, sig, tries = 2)
    nulls[b] <- 2*(g1$val - g2$val)
  }
  p_blrt <- (1 + sum(nulls >= lrt)) / (B + 1)
  c(bic_fires = as.numeric(bic2 < bic1), blrt_p = p_blrt, lrt = lrt)
}

run_cell <- function(n, two_class, ratio, ndraw, B) {
  sep <- ratio * 0.3 * 12.93 / log(10)
  bic <- 0; blrt <- 0
  for (i in seq_len(ndraw)) {
    sig <- het(n)
    if (two_class) { n2 <- rbinom(1, n, F2S); mu <- c(rep(0, n - n2), rep(sep, n2)) }
    else           { mu <- rep(0, n) }
    x <- rnorm(n, mu, sqrt(SINT^2 + sig^2))
    d <- decide(x, sig, B)
    bic <- bic + as.numeric(d["bic_fires"]); if (as.numeric(d["blrt_p"]) <= 0.05) blrt <- blrt + 1
  }
  c(bic = bic/ndraw, blrt = blrt/ndraw)
}

if (length(args) >= 1) {
  # ---- chunk mode, for shell-level parallelism over many draws
  two <- (MODE == "power")
  for (i in seq_len(NDRAW)) {
    sig <- het(N)
    if (two) { n2 <- rbinom(1, N, F2S); mu <- c(rep(0, N - n2), excess(n2)) }
    else     { mu <- rep(0, N) }
    x <- rnorm(N, mu, sqrt(SINT^2 + sig^2))
    d <- decide(x, sig, B)
    cat(sprintf("%d\t%s\t%.2f\t%.1f\t%.0f\t%.4f\t%.3f\n", N, MODE, RATIO, CONC,
                d["bic_fires"], d["blrt_p"], d["lrt"]))
  }
} else {
  # ---- standalone: a reduced grid that reproduces the qualitative result in a few minutes
  NB <- 39; ND <- 8; NN <- 1000
  cat("=== reduced grid: n = 1000, ", ND, " draws per cell, B = ", NB, " bootstraps\n\n", sep="")
  cat("   case                       BIC fires   BLRT fires\n")
  r <- run_cell(NN, FALSE, 0.5, ND, NB)
  cat(sprintf("   one class (null)         %10.3f %12.3f\n", r["bic"], r["blrt"]))
  r <- run_cell(NN, TRUE,  0.5, ND, NB)
  cat(sprintf("   far side at half         %10.3f %12.3f\n", r["bic"], r["blrt"]))
  r <- run_cell(NN, TRUE,  0.2, ND, NB)
  cat(sprintf("   far side at a fifth      %10.3f %12.3f\n", r["bic"], r["blrt"]))

  cat("\n=== the full grid, 465 draws, B = 99, run 2026-09-21 ===\n\n")
  cat("   Reproduce with, for each cell:\n")
  cat("     Rscript blrt_vs_bic.R <n> <power|null> <seed> <ndraws> 99 <ratio>\n\n")
  cat("      far side     n      draws   BIC fires   BLRT fires\n")
  cat("      (null)       300       70       0.000        0.057\n")
  cat("      (null)      1000       70       0.000        0.057\n")
  cat("      (null)      3000       50       0.000        0.040\n")
  cat("      half         300       60       0.000        0.250\n")
  cat("      half        1000       60       0.033        0.567\n")
  cat("      half        3000       45       0.356        0.978\n")
  cat("      a fifth      300       40       0.000        0.050\n")
  cat("      a fifth     1000       40       0.000        0.075\n")
  cat("      a fifth     3000       30       0.000        0.033\n")

  cat("\n=== flatly\n\n")
  cat("  BLRT IS VALID HERE. On one-class data it fires 0.057, 0.057 and 0.040 against a\n")
  cat("  nominal 0.05, so the power below is not bought by being anticonservative. That had\n")
  cat("  to be checked first: the naive estimator deconvolved_mixture.R killed was killed by\n")
  cat("  its null, not by its power.\n\n")
  cat("  BLRT BEATS BIC BY A LOT at the separation 4.3 quotes. 0.978 against 0.356 at three\n")
  cat("  thousand objects, and 0.567 against 0.033 at one thousand. BIC is the wrong tool for\n")
  cat("  choosing between one component and two: the mixture weight sits on the boundary of\n")
  cat("  the parameter space, the usual asymptotics do not hold, and bootstrapping the null\n")
  cat("  distribution of the LRT is the standard fix in latent-class work. That is the\n")
  cat("  borrowable calculation, and it is borrowed from psychometrics rather than invented.\n\n")
  cat("  THE FIFTH-RATE NEGATIVE SURVIVES, and this is the useful part. 4.3 said the test\n")
  cat("  fails at every size tried when the far side supplies a fifth, and that claim came\n")
  cat("  from the NAIVE analysis the deconvolution superseded. Re-run with the valid and\n")
  cat("  more powerful estimator it still fails: 0.050, 0.075, 0.033 at n = 300, 1000, 3000,\n")
  cat("  which is the nominal false-positive rate and nothing more. A 0.337 dex separation\n")
  cat("  inside 0.40 dex of scatter is not recoverable by a better test, so the negative is\n")
  cat("  about the geometry and not about the estimator.\n\n")
  cat("  WHAT THIS CHANGES IN 4.3: the sample size. The test is feasible at one to three\n")
  cat("  thousand objects with the right estimator, not 'well beyond three thousand'.\n\n")
  cat("=== the duty cycle, put back where the physics has it (CONC = 4)\n\n")
  cat("  Everything above holds every class-2 object at the MEAN duty cycle, which is what\n")
  cat("  deconvolved_mixture.R did. The physics is per-object: the excess is exp((k-1) lam N1),\n")
  cat("  so with lam ~ Beta(mean 0.3, conc 4) the excess has mean 0.844 dex and sd 0.577, and\n")
  cat("  27 per cent of class 2 sits below one intrinsic scatter. I expected that to cost\n")
  cat("  power. It does the opposite.\n\n")
  cat("      n      BIC fixed -> spread     BLRT fixed -> spread\n")
  cat("      300      0.000 -> 0.575          0.250 -> 0.925\n")
  cat("     1000      0.033 -> 1.000          0.567 -> 1.000\n")
  cat("     3000      0.356 -> 1.000          0.978 -> 1.000\n\n")
  cat("  The null is untouched, because CONC never enters the one-class branch; pooled over\n")
  cat("  90 draws it is 0.067 at n=300 and 0.078 at n=1000, both consistent with 0.05.\n\n")
  cat("  WHY IT HELPS. A spread of duty cycles gives class 2 a large VARIANCE as well as a\n")
  cat("  shifted mean, and a two-component mixture picks up a broad component beside a narrow\n")
  cat("  one more readily than a small mean shift. So the fixed-duty-cycle numbers quoted in\n")
  cat("  4.3 are CONSERVATIVE, and the test is within reach of a few hundred objects rather\n")
  cat("  than a few thousand.\n\n")
  cat("  AND THE CAVEAT THAT COMES WITH IT, which is why 4.3 keeps the conservative number as\n")
  cat("  its headline. What the test then detects is partly the BREADTH of the second class and\n")
  cat("  not a clean mean shift, and breadth has other astrophysical causes: a subpopulation\n")
  cat("  with larger intrinsic scatter, or with worse measurements, would also fire it. A mean\n")
  cat("  shift is more specific evidence for two classes than a variance difference is.\n")
}

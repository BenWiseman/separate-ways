# two_sided_growth_exponent.R ends by proposing bimodality in the M_BH-M_star residual
# as a duty-cycle-free test of the two-sided class. That was asserted, not shown, and it
# has an obvious way to fail: the two-sided sequence carries the factor exp((k-1) lam N1),
# so objects with a low duty cycle sit at ZERO excess, on top of the one-sided sequence.
# If the duty-cycle distribution has weight near zero the two sequences merge and there
# is no bimodality to find. Test it before it goes in a paper.
#
# Base R only, per the repo rule. The mixture fit is a hand-rolled EM.

set.seed(20260921)
N1 <- 12.93          # one-sided e-folds, z=20 to z=7 (two_sided_growth_exponent.R)

# --- a two-component Gaussian mixture by EM, and a one-component fit, compared by BIC
fit1 <- function(x) {
  mu <- mean(x); s <- sd(x)
  ll <- sum(dnorm(x, mu, s, log=TRUE))
  list(ll=ll, k=2, bic=-2*ll + 2*log(length(x)))
}
fit2 <- function(x, tries=6) {
  best <- NULL
  for (t in seq_len(tries)) {
    q <- quantile(x, c(0.25, 0.75)) + rnorm(2, 0, 0.1*sd(x))
    mu <- as.numeric(q); s <- rep(sd(x)/2, 2); p <- 0.5
    for (it in 1:300) {
      d1 <- p*dnorm(x, mu[1], s[1]); d2 <- (1-p)*dnorm(x, mu[2], s[2])
      tot <- d1 + d2; tot[tot <= 0] <- 1e-300
      w <- d1/tot
      p  <- mean(w)
      if (p < 1e-4 || p > 1-1e-4) break
      mu <- c(sum(w*x)/sum(w), sum((1-w)*x)/sum(1-w))
      s  <- c(sqrt(sum(w*(x-mu[1])^2)/sum(w)), sqrt(sum((1-w)*(x-mu[2])^2)/sum(1-w)))
      s[s < 1e-3] <- 1e-3
    }
    ll <- sum(log(pmax(p*dnorm(x,mu[1],s[1]) + (1-p)*dnorm(x,mu[2],s[2]), 1e-300)))
    if (is.finite(ll) && (is.null(best) || ll > best$ll))
      best <- list(ll=ll, mu=mu, s=s, p=p)
  }
  best$bic <- -2*best$ll + 5*log(length(x))
  best
}

# --- draw a population
draw <- function(n, f2s, kappa, lam_mean, lam_conc, sig_int) {
  n2 <- rbinom(1, n, f2s)
  a <- lam_mean*lam_conc; b <- (1-lam_mean)*lam_conc
  lam <- rbeta(n2, a, b)
  excess2 <- (kappa-1)*lam*N1/log(10)          # dex
  c(rnorm(n-n2, 0, sig_int), rnorm(n2, excess2, sig_int))
}

cat("=== 1. does the two-sided sequence separate at all?\n\n")
cat("   Mean excess in dex for the two-sided objects, before scatter:\n\n")
cat("      kappa   <lambda>   mean excess (dex)   as a factor\n")
for (k in c(1.2, 1.5, 2.0)) for (lm in c(0.1, 0.3, 0.6)) {
  e <- (k-1)*lm*N1/log(10)
  cat(sprintf("   %7.2f %10.2f %19.2f %13.3g\n", k, lm, e, 10^e))
}
cat("\n   Against an intrinsic M_BH-M_star scatter of 0.3-0.5 dex, separations below\n")
cat("   about 0.5 dex are hopeless and above about 1.5 dex are obvious.\n")

cat("\n=== 2. detection rate, by simulation\n\n")
cat("   Fraction of 200 mock catalogues in which a two-component mixture beats a\n")
cat("   single Gaussian on BIC. sigma_int = 0.4 dex, f_2s = 0.3, Beta duty cycle.\n\n")
sig_int <- 0.4; f2s <- 0.3; conc <- 4
cat("      kappa  <lambda>   mean sep (dex)      n=100   n=300  n=1000\n")
for (k in c(1.2, 1.5, 2.0)) for (lm in c(0.3, 0.6)) {
  sep <- (k-1)*lm*N1/log(10)
  cat(sprintf("   %7.2f %9.2f %16.2f", k, lm, sep))
  for (n in c(100, 300, 1000)) {
    hit <- 0
    for (r in 1:200) {
      x <- draw(n, f2s, k, lm, conc, sig_int)
      if (fit2(x)$bic < fit1(x)$bic) hit <- hit + 1
    }
    cat(sprintf(" %7.2f", hit/200))
  }
  cat("\n")
}

cat("\n=== 3. the null: how often does BIC cry bimodal when there is none?\n\n")
cat("   Same test on a pure one-sided population (f_2s = 0), which must not fire.\n\n")
cat("        n      false-positive rate\n")
for (n in c(100, 300, 1000)) {
  hit <- 0
  for (r in 1:200) {
    x <- rnorm(n, 0, sig_int)
    if (fit2(x)$bic < fit1(x)$bic) hit <- hit + 1
  }
  cat(sprintf("   %6d %22.3f\n", n, hit/200))
}

cat("\n=== 3b. the check that matters: measurement error, not just intrinsic scatter\n\n")
cat("   Section 2 used sigma = 0.4 dex, which is the INTRINSIC scatter of the relation.\n")
cat("   The observable is a residual built from two measured quantities, and both are\n")
cat("   badly measured in these objects: virial M_BH from broad lines carries ~0.5 dex of\n")
cat("   systematic, SED-fitted M_star ~0.3 dex. Those add in quadrature ON TOP of the\n")
cat("   intrinsic scatter, and the paper now quotes a sample size, so check it.\n\n")
sig_tot <- function(si, eb, es) sqrt(si^2 + eb^2 + es^2)
cat("      sigma_int   err(M_BH)  err(M_star)   total sigma\n")
for (v in list(c(0.4,0,0), c(0.4,0.3,0.2), c(0.4,0.5,0.3))) {
  cat(sprintf("   %10.2f %11.2f %12.2f %13.3f\n", v[1], v[2], v[3], sig_tot(v[1],v[2],v[3])))
}
cat("\n   Detection rate at kappa=1.5, <lambda>=0.3 (separation 0.84 dex), f_2s=0.3,\n")
cat("   100 draws per cell, against total sigma rather than intrinsic:\n\n")
cat("      total sigma    sep/sigma    n=300   n=1000   n=3000\n")
for (v in list(c(0.4,0,0), c(0.4,0.3,0.2), c(0.4,0.5,0.3))) {
  st <- sig_tot(v[1],v[2],v[3]); sep <- 0.5*0.3*N1/log(10)
  cat(sprintf("   %12.3f %12.2f", st, sep/st))
  for (n in c(300, 1000, 3000)) {
    hit <- 0
    for (r in 1:100) {
      x <- draw(n, 0.3, 1.5, 0.3, 4, st)
      if (fit2(x, tries=4)$bic < fit1(x)$bic) hit <- hit + 1
    }
    cat(sprintf(" %8.2f", hit/100))
  }
  cat("\n")
}

cat("\n=== 4. flatly\n\n")
cat("  The proposal survives, conditionally, and the condition is the thing to quote.\n")
cat("  A two-sided fraction of a third with kappa at or above 1.5 and a duty cycle at\n")
cat("  or above 0.3 is detectable in a few hundred objects ONLY if the residual carries\n")
cat("  the relation's intrinsic scatter alone. Section 3b prices the real case: fold in\n")
cat("  the measurement errors on a virial black hole mass and a fitted stellar mass and\n")
cat("  the power at n=300 falls from 1.00 to 0.16, needing a few thousand objects\n")
cat("  instead. The wall is measurement error, not sample size, and that is a factor of\n")
cat("  ten. At kappa=1.2, which is\n")
cat("  what a face-value reading of the overmassive ratios gives, the separation is\n")
cat("  under half a dex and the test does not fire at any sample size reached here:\n")
cat("  the two sequences are inside the intrinsic scatter and no amount of data\n")
cat("  separates them, because the confusion is systematic and not statistical.\n\n")
cat("  So the honest statement for the paper is a conditional one. Bimodality is a\n")
cat("  real, duty-cycle-robust signature of the two-sided class, and it is reachable\n")
cat("  only in the half of parameter space where the far side contributes at least\n")
cat("  half of what ours does. That is worth writing down because it is the first\n")
cat("  test of the two-sided class that does not need a ringdown.\n\n")
cat("  NEXT ROUTE, two of them, and the second is the one worth asking an expert.\n")
cat("  (a) Reduce the scatter rather than adding objects. The fold's excess is in M_BH\n")
cat("      alone and M_BH-sigma is tighter than M_BH-M_star, so running the same test\n")
cat("      against a velocity dispersion buys the difference in the two relations'\n")
cat("      scatters. Whether that is enough is a catalogue question.\n")
cat("  (b) Deconvolve instead of inflating. Section 3b folds measurement error into the\n")
cat("      component widths, which is the naive thing to do and throws away the fact\n")
cat("      that the per-object errors are KNOWN. A mixture fitted with per-object error\n")
cat("      variances -- measurement-error mixture modelling, standard in psychometrics\n")
cat("      and behavioural phenotyping -- recovers some of the factor of ten, and how\n")
cat("      much is an open question this script cannot answer. That is a well-posed ask\n")
cat("      for a quantitative methodologist and it has a number at the end of it.\n")

# How many events does the falsifier ACTUALLY need?
#
# calibrated_falsifier.R tests the endpoint by counting events above a calibrated threshold and
# firing if any exceed it. That is honest and assumption-light, and it throws away every event
# below the threshold, which at a hundred events is almost all of them. The endpoint is a
# boundary-of-support problem, and for a simple alternative Neyman-Pearson says the likelihood
# ratio on the WHOLE sample is the most powerful test at the same size. Nobody had asked how much
# is being left on the table, so this asks.
#
# Same setup as calibrated_falsifier.R, unchanged: E^-2 above 50 PeV truncated at the endpoint,
# lognormal response of width 0.3 in ln E, model endpoint 245.8 PeV, alternatives at 2x, 3x, 5x.
# Both tests are calibrated by Monte Carlo to the SAME 5 per cent sample-wide false alarm, so the
# comparison is like for like.
#
# Base R only, per the repo rule: no package is loaded.

set.seed(20260921)
E_m <- 245.8; Ecut <- 50; sg <- 0.3

# ---- reconstructed-energy density for a source population truncated at Emax.
# The source is E^-2 dE on [Ecut, Emax]; on a log grid that integrand is exp(-lnE). The response
# is lognormal in E, so the reconstructed density is a Gaussian convolution in ln E.
make_logdens <- function(Emax, n_t = 1400, n_o = 1800) {
  lt <- seq(log(Ecut), log(Emax), length.out = n_t)
  tw <- c(diff(lt)[1]/2, rep(diff(lt)[1], n_t-2), diff(lt)[1]/2)   # trapezoid in ln E
  src <- exp(-lt)                                                   # E^-2 * E
  Z <- sum(tw * src)
  lo <- seq(log(Ecut) - 5*sg, log(Emax) + 6*sg, length.out = n_o)
  # p(lo) as a density in ln E_obs, then convert to a density in E_obs by dividing by E_obs
  dens_ln <- sapply(lo, function(L)
    sum(tw * src * dnorm(L - lt, 0, sg)) / Z)
  list(lo = lo, ln_p = log(pmax(dens_ln, 1e-300)),
       check = sum(c(diff(lo)[1]/2, rep(diff(lo)[1], n_o-2), diff(lo)[1]/2) * dens_ln))
}
logp_fun <- function(dd) function(E) {
  approx(dd$lo, dd$ln_p, xout = log(E), rule = 2)$y - log(E)   # density in E_obs
}

# ---- draw a sample: E^-2 on [Ecut, Emax] by inverse CDF, then smear
draw <- function(n, Emax) {
  u <- runif(n)
  Et <- 1/(1/Ecut - u*(1/Ecut - 1/Emax))
  Et * exp(rnorm(n, 0, sg))
}

# ---- validate the densities before using them
D <- list()
for (r in c(1, 2, 3, 5)) D[[as.character(r)]] <- make_logdens(r*E_m)
for (r in names(D))
  stopifnot(abs(D[[r]]$check - 1) < 2e-3)        # normalised in ln E_obs
cat(sprintf("density normalisation, r = 1,2,3,5: %s\n",
            paste(sprintf("%.5f", sapply(D, function(d) d$check)), collapse=", ")))

# the drawn sample's mean log-likelihood must be highest under its OWN hypothesis
set.seed(7)
for (r in c(2, 3, 5)) {
  x <- draw(20000, r*E_m)
  ll_true <- mean(logp_fun(D[[as.character(r)]])(x))
  ll_null <- mean(logp_fun(D[["1"]])(x))
  stopifnot(ll_true > ll_null)
}
cat("likelihood check: each sample scores highest under the hypothesis that generated it\n\n")

# ---- the two tests, both calibrated by Monte Carlo to 5 per cent sample-wide
NMC <- 4000; ALPHA <- 0.05
lp0 <- logp_fun(D[["1"]])

# threshold test, as calibrated_falsifier.R builds it
p_above <- function(thr, Emax) {
  f <- function(E) (1/E^2) * pnorm((log(E/thr))/sg)
  integrate(f, Ecut, Emax, subdivisions=4000)$value /
    integrate(function(E) 1/E^2, Ecut, Emax)$value
}
thr_for <- function(N, alpha=ALPHA)
  uniroot(function(t) p_above(t, E_m) - (1-(1-alpha)^(1/N)), c(E_m, 60*E_m), tol=1e-10)$root

cat("   both tests at a 5 per cent sample-wide false alarm, model endpoint 245.8 PeV\n\n")
cat("                        threshold test              likelihood ratio\n")
cat("      N       r=2    r=3    r=5        r=2    r=3    r=5\n")
for (N in c(10, 30, 100, 300)) {
  th <- thr_for(N)
  pw_thr <- sapply(c(2,3,5), function(r) 1-(1-p_above(th, r*E_m))^N)

  pw_lrt <- numeric(3)
  for (k in seq_along(c(2,3,5))) {
    r <- c(2,3,5)[k]; lp1 <- logp_fun(D[[as.character(r)]])
    null <- replicate(NMC, { x <- draw(N, E_m);     sum(lp1(x) - lp0(x)) })
    alt  <- replicate(NMC, { x <- draw(N, r*E_m);   sum(lp1(x) - lp0(x)) })
    crit <- quantile(null, 1-ALPHA, names = FALSE)
    pw_lrt[k] <- mean(alt > crit)
  }
  cat(sprintf("   %5d %8.3f %6.3f %6.3f %10.3f %6.3f %6.3f\n",
              N, pw_thr[1], pw_thr[2], pw_thr[3], pw_lrt[1], pw_lrt[2], pw_lrt[3]))
}

# ---- where does the likelihood ratio reach the threshold test's power at N=100?
th100 <- thr_for(100); target <- 1-(1-p_above(th100, 2*E_m))^100
lp1 <- logp_fun(D[["2"]])
cat(sprintf("\n   the threshold test reaches %.3f against a doubled endpoint at N = 100.\n", target))
cat("   the likelihood ratio, finer scan in N:\n\n      N     power\n")
cross <- NA
for (N in c(30, 40, 45, 50, 60)) {
  null <- replicate(NMC, { x <- draw(N, E_m);     sum(lp1(x) - lp0(x)) })
  alt  <- replicate(NMC, { x <- draw(N, 2*E_m);   sum(lp1(x) - lp0(x)) })
  pw <- mean(alt > quantile(null, 1-ALPHA, names = FALSE))
  cat(sprintf("   %4d %9.3f%s\n", N, pw, if (is.na(cross) && pw >= target) "   <- crosses here" else ""))
  if (is.na(cross) && pw >= target) cross <- N
}
cat(sprintf("\n   so the whole-sample likelihood buys a factor of about %.1f in events.\n",
            100/cross))

cat("\n=== flatly\n\n")
cat("  Read the two halves against each other. The likelihood ratio uses every event; the\n")
cat("  threshold test uses only those above a cut that, at a hundred events, almost nothing\n")
cat("  reaches. If the right-hand block reaches at ten or thirty events the power the left\n")
cat("  block needs a hundred for, the falsifier is cheaper than 3.2 says it is, and by how\n")
cat("  much is the number to quote.\n\n")
cat("  WHAT IT COSTS. The likelihood ratio assumes the E^-2 shape and the response width are\n")
cat("  known; the threshold test does not. That is a real trade and 3.2 should say so rather\n")
cat("  than quietly taking the better number. A spectrum that is not E^-2, or a response tail\n")
cat("  fatter than lognormal, degrades the likelihood ratio and leaves the threshold test\n")
cat("  alone. The honest reading is that the threshold test is the conservative falsifier and\n")
cat("  the likelihood ratio is what is achievable once the instrument is understood.\n")

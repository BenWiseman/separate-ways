#!/usr/bin/env Rscript
# ringdown_null_prediction.R
#
# The echo idea is dead and the kill is the paper's own. The pinch at Kruskal
# V = 59.21 sits INSIDE the horizon, and A.14 already establishes that anything
# there "lies in the causal future of both exteriors and therefore reaches
# neither". Separately the adopted implementation sets the seam weight to zero,
# so there is no horizon reflectivity and nothing to echo from. Two independent
# reasons, both already in the manuscript. Do not compute an exterior observable
# from a causally sealed region; that error was made once already in this project
# with an MeV gamma-ray line.
#
# WHAT SURVIVES IS THE NULL ITSELF. Every gravitational-wave signature the fold
# could produce (horizon reflectivity, the l -> l+/-1 multipole mixing of A.11,
# the omega -> omega* spectral action of A.12) requires a NONZERO seam weight,
# and the implementation this paper defends sets it to zero. So the paper makes a
# sharp prediction: ordinary Kerr ringdown, no echo, no damping-time shift. That
# is falsifiable, and there is an instrument taking data against it now.

cat("=== 1. the prediction ===\n")
cat("  transparent seam (adopted): delta_tau_220 = 0 EXACTLY, no echo, pure Kerr\n")
cat("  reflecting seam (worked out separately): nonzero, spin-dependent\n")

cat("\n=== 2. what the existing measurement says ===\n")
# GWTC-3, as quoted in the companion: 90% credible intervals
est <- list(
  list(name="posterior multiplication", mu=0.14, lo=0.11, hi=0.11),
  list(name="hierarchical combination",  mu=0.13, lo=0.22, hi=0.21))
for (e in est) {
  # 90% credible -> approximate 1-sigma by dividing by 1.645
  sig <- mean(c(e$lo, e$hi))/1.645
  z   <- e$mu/sig
  cat(sprintf("  %-26s delta_tau = %.2f (+%.2f/-%.2f, 90%%) -> 1-sigma %.3f, null at %.1f sigma\n",
              e$name, e$mu, e$hi, e$lo, sig, z))
}

cat("\n=== 3. validation: the check must be able to fail ===\n")
z1 <- 0.14/(0.11/1.645); z2 <- 0.13/(0.215/1.645)
stopifnot(z1 > z2)            # the tighter interval must give the larger tension
cat(sprintf("  tighter interval gives larger tension (%.1f vs %.1f sigma), as it must\n", z1, z2))
bad <- tryCatch({ stopifnot(z1 < z2); TRUE }, error=function(e) FALSE)
cat(sprintf("  reversed assertion fails: %s\n", if(!bad) "yes" else "NO - BLIND"))

cat("\n=== 4. how much better does a measurement have to get? ===\n")
cat("  If the true deviation were the central value 0.14, the exposure needed to\n")
cat("  reach a given significance against the transparent seam's exact zero:\n\n")
sig_now <- 0.11/1.645
cat(sprintf("  %10s %14s %16s\n", "target", "sigma needed", "events vs now"))
for (target in c(3, 5)) {
  need <- 0.14/target
  scale <- (sig_now/need)^2
  cat(sprintf("  %10s %14.3f %16.1f x\n", paste0(target, " sigma"), need, scale))
}

cat(sprintf("
=== 5. flatly ===

  The echo route is dead and was dead before it was tried: the pinch is inside a
  horizon and the paper's own A.14 seals it, while the adopted seam reflects
  nothing. No exterior observable follows from the 59.21 depth and none should be
  claimed.

  What the paper does have is a null prediction with a live instrument. Every
  gravitational-wave signature the fold could produce needs a nonzero seam weight,
  and this implementation sets it to zero, so it predicts an ordinary Kerr
  ringdown: no echo, and a damping-time deviation of exactly zero.

  The current measurement is the interesting part. GWTC-3 gives delta_tau_220 =
  0.14 +/- 0.11 by posterior multiplication, which is %.1f sigma from zero, and
  0.13 +/- 0.22 hierarchically, which is %.1f sigma and consistent. The LVK
  analysis itself does not claim a violation and discusses positive bias. So the
  transparent seam is not excluded, and it is not comfortable either.

  That makes this a rare thing for the paper: a place where the model says
  something sharp about data being taken right now. If the damping-time excess
  firms up, the implementation defended here is excluded and the reflecting branch
  the paper works out separately takes over, with its own spin-dependent
  signature. If the ringdown settles onto pure Kerr, this implementation stands
  and the reflecting alternatives are the ones in trouble.

  Reaching three sigma on the current central value needs about %.0f times the
  present exposure, which is an O4/O5 question rather than a next-decade one.
", z1, z2, (sig_now/(0.14/3))^2))

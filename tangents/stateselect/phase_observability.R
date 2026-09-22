# F1 asks what fixes the residual phase mu(p). Four selection principles have died trying.
# Inference asks a question F1 does not: how much could the ABUNDANCE ever say about mu?
#
# The family is n_mu(x) = (1 - sqrt(1-P) cos(mu - mu_*))/2 with P = exp(-x^2), and the observable
# is I(mu) = (1/pi^2) int x^2 n_mu dx, which fixes the mass. Differentiate under the integral.
P     <- function(x) exp(-x^2)
D     <- function(x) sqrt(1-P(x))                     # band half-width, n_max - n_min = D
w     <- function(f,lo=0,hi=30) integrate(function(x) x^2*f(x), lo, hi, subdivisions=8000)$value
# The phase excursion must DECAY, or the state is inadmissible: a constant offset d sends
# n -> sin^2(d/2) at large x, a floor that makes int x^2 n dx diverge. A first version of this
# file used a constant offset and produced nonsense (I/I_0 = 45 at d = 0.05). Use a decaying
# profile, which is what A.18's ultraviolet condition actually permits.
dprof <- function(x, a) a*exp(-x^2)
I_of  <- function(a) w(function(x) (1 - D(x)*cos(dprof(x,a)))/2)/pi^2
I0    <- I_of(0)
cat(sprintf("  I(mu_*) = %.10f   (5.1 uses 0.0127597)\n\n", I0))

cat("=== 1. the abundance is stationary at the minimum, so it is blind to first order ===\n")
for (d in c(1e-4, 1e-3, 1e-2)) {
  num <- (I_of(d) - I_of(-d))/(2*d)                   # central first derivative
  cat(sprintf("   step %.0e : dI/d(mu) = %+.3e\n", d, num))
}
cat("   The first derivative vanishes: mu_* is a stationary point of the observable, not just\n")
cat("   of the occupation. Any abundance measurement is insensitive to the phase at leading order.\n")

cat("\n=== 2. so the leading sensitivity is quadratic. How big is it? ===\n")
Ipp <- (I_of(1e-3) - 2*I_of(0) + I_of(-1e-3))/1e-6    # I'' at a=0, numerically
cat(sprintf("   I''(mu_*) = %.10f, so dI/I = (I''/2I0) d^2 = %.6f d^2\n", Ipp, Ipp/(2*I0)))
cat("\n     abundance precision   smallest detectable phase excursion d (rad)\n")
for (prec in c(1e-2, 1e-3, 1e-4, 1e-5)) {
  d <- sqrt(prec*2*I0/Ipp)
  cat(sprintf("   %19.0e %36.4f\n", prec, d))
}
cat("\n   Omega_DM h^2 is known to about one per cent, which leaves the phase unconstrained\n")
cat("   until |mu - mu_*| exceeds roughly the value in the first row.\n")

cat("\n=== 3. what that does to the mass, which is the quantity anyone quotes ===\n")
M0 <- 491.6
cat("   amplitude a of a exp(-x^2)   I/I_0      mass (PeV)      shift\n")
for (d in c(0.1, 0.3, 1.0, 2.0, 3.0)) {
  r <- I_of(d)/I0
  cat(sprintf("   %18.2f %11.6f %14.2f %11.2f\n", d, r, M0*r^(-2/5), M0*r^(-2/5)-M0))
}
cat("\n  FLATLY. This does not fix the phase and does not pretend to. What it says is that the\n")
cat("  abundance cannot fix it either, and by how much: the observable is quadratically flat at\n")
cat("  the adopted point, so a one-per-cent abundance leaves a whole band of phases\n")
cat("  indistinguishable. F1 is not merely unsolved, it is unsolvable BY THIS OBSERVABLE.\n")
cat("\n  NEXT ROUTE. A quantity linear in the phase would beat a quadratic one. R, the decoherence\n")
cat("  normalisation, is built from the min-entropy rather than the mean, so its derivative at\n")
cat("  mu_* need not vanish. Check whether dR/dmu is non-zero there; if it is, the decoherence\n")
cat("  clock constrains the phase where the abundance cannot, and 4.2 becomes the sharper probe.\n")

cat("\n=== 4. the coherence worth noticing ===\n")
a_inv <- uniroot(function(a) I_of(a)/I0 - 1.01, c(1e-4, 2))$root
cat(sprintf("   phase amplitude invisible at one per cent abundance precision: a = %.4f rad\n", a_inv))
cat(sprintf("   mass shift it hides: %.2f PeV, against a quoted width of 2.0 PeV\n",
    491.6 - 491.6*(I_of(a_inv)/I0)^(-2/5)))
cat("   Not a coincidence, since both trace to the same one per cent, but useful: the quoted\n")
cat("   width already covers the state freedom the abundance cannot resolve.\n")

cat("\n=== 5. the next route was tested and it fails, for a reason that generalises ===\n")
Rn <- function(a) {
  num <- w(function(x){u <- pmin(pmax((1-D(x)*cos(dprof(x,a)))/2, 1e-14), 1-1e-14)
                        -log(pmax(u, 1-u))})
  num/w(function(x) (1-D(x)*cos(dprof(x,a)))/2)
}
cat("   R is built from the min-entropy rather than the mean, so it might have been linear.\n")
for (h in c(1e-3, 1e-2, 5e-2))
  cat(sprintf("     dR/da at a = 0, step %.0e : %+.4e\n", h, (Rn(h)-Rn(-h))/(2*h)))
cat("\n   It is not. And the reason is general rather than accidental: mu_* minimises the\n")
cat("   occupation POINTWISE, so dn/da vanishes at every x when a = 0,\n")
cat("     dn/da = D(x) sin(a e^{-x^2}) e^{-x^2} / 2  ->  0.\n")
cat("   Any differentiable functional F[n] therefore has dF/da = int (dF/dn)(dn/da) = 0 there.\n")
cat("\n  SO: no observable built from the occupation can constrain the residual phase at first\n")
cat("  order. Not the abundance, not the decoherence clock, not any other functional of n.\n")
cat("  F1 is not merely unsolved; the class of measurements the paper has access to cannot\n")
cat("  solve it. Constraining the phase needs something sensitive to it directly, an\n")
cat("  interference observable rather than an occupation, and the cosmology offers none.\n")

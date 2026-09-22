# n_* = (1-sqrt(1-P))/2 and C = sqrt(P), so n_* = (1-sqrt(1-C^2))/2: a function of C alone.
# The asymptotic Fermi-Dirac form with mu = log 4 is n_FD = 1/(4/P + 1) = C^2/(4+C^2):
# ALSO a function of C alone. So the departure from thermality is an exact function of the
# entanglement. Get its small-C behaviour in closed form and check it.
nstar <- function(C) (1 - sqrt(1-C^2))/2
nFD   <- function(C) C^2/(4+C^2)
dev   <- function(C) nstar(C) - nFD(C)

cat("       C        n_*(C)        n_FD(C)         deviation      deviation/C^4\n")
for (C in c(1,0.9,0.5,0.2,0.1,0.05,0.01,1e-3))
  cat(sprintf("   %8.4f %13.9f %13.9f %15.6e %15.9f\n", C, nstar(C), nFD(C), dev(C), dev(C)/C^4))
# below C ~ 1e-3 the difference of two nearly equal doubles loses all its digits; the
# ratio there is floating-point noise, not a change in the power law, so the table stops.

cat("\n  Series: n_* = C^2/4 + C^4/16 + O(C^6);  n_FD = C^2/4 - C^4/16 + O(C^6).\n")
cat("  So the leading terms agree and the deviation is C^4/8 + O(C^6). Check the coefficient:\n")
for (C in c(0.1,0.05,0.01,1e-3)) cat(sprintf("     C = %.4f : deviation/C^4 = %.10f   (1/8 = 0.125)\n", C, dev(C)/C^4))

cat("\n  And at maximal entanglement:\n")
cat(sprintf("     C = 1 : n_* = %.6f, n_FD = %.6f, deviation = %.6f\n", nstar(1), nFD(1), dev(1)))
cat(sprintf("     the relic is %.1f per cent more occupied than a thermal state of the same\n",
            100*dev(1)/nFD(1)))
cat("     fugacity would be, exactly where the pair is a Bell state.\n")

# where does the production integral actually live? thermal region or entangled region?
P <- function(x) exp(-x^2); Cx <- function(x) sqrt(P(x))
w  <- function(x) x^2*nstar(Cx(x))
tot <- integrate(w, 0, 40, rel.tol=1e-12)$value
for (Cc in c(0.5, 0.1, 0.01)) {
  xc <- sqrt(-2*log(Cc))
  fr <- integrate(w, 0, xc, rel.tol=1e-12)$value/tot
  cat(sprintf("\n  fraction of Int x^2 n dx from modes with C > %.2f (x < %.3f): %.4f", Cc, xc, fr))
}
cat("\n\n  FLATLY: the state is NOT a thermofield double. It is thermal only where it carries no\n")
cat("  entanglement, and the departure is an exact function of the concurrence, C^4/8 to\n")
cat("  leading order. Thermality and entanglement are traded off mode by mode.\n")
# Is the two-sheet state a thermofield double? If the fold is the modular conjugation and the
# global state is its canonical purification, tracing out the partner sheet must leave OUR
# sheet in a Gibbs state, i.e. the occupation must be Fermi-Dirac:  n = 1/(e^{beta w}+1),
# equivalently  log((1-n)/n) = beta*w  must be LINEAR in whatever the energy is.
# Test it rather than assert it.
P  <- function(x) exp(-x^2)
ns <- function(x) (1 - sqrt(1-P(x)))/2
lg <- function(x) log((1-ns(x))/ns(x))     # = beta*w if the state is thermal

x <- seq(0.05, 6, by=0.01)
y <- lg(x)
cat("  If thermal in a RELATIVISTIC energy w ∝ x, log((1-n)/n) must be linear in x.\n")
f1 <- lm(y ~ x);  cat(sprintf("     fit vs x   : R^2 = %.6f   max resid = %.4f\n",
                              summary(f1)$r.squared, max(abs(resid(f1)))))
cat("  If thermal in a NON-RELATIVISTIC energy w ∝ x^2, it must be linear in x^2.\n")
x2 <- x^2; f2 <- lm(y ~ x2); cat(sprintf("     fit vs x^2 : R^2 = %.8f   max resid = %.4e\n",
                              summary(f2)$r.squared, max(abs(resid(f2)))))
cat(sprintf("     slope = %.10f   intercept = %.10f   (log 4 = %.10f)\n",
            coef(f2)[2], coef(f2)[1], log(4)))

cat("\n  Exact statement, checked mode by mode against  n = 1/(e^{x^2 + mu}+1):\n")
cat("        x        n_*(x)        Fermi-Dirac(x^2, mu=log4)      difference\n")
for (xx in c(0.1,0.5,1,1.5,2,3,4,5)) {
  fd <- 1/(exp(xx^2)*4 + 1)
  cat(sprintf("   %6.2f  %14.10f  %22.10f  %14.3e\n", xx, ns(xx), fd, ns(xx)-fd))
}
cat("\n  And the exact algebra: n_* = (1-sqrt(1-e^{-x^2}))/2, so\n")
cat("     log((1-n_*)/n_*) - x^2  ->  log 4 only as x -> infinity.\n")
# stop at x = 5: beyond it n_* underflows to zero in double precision and the ratio reports
# Inf, which is the floating-point floor and not a statement about the state.
for (xx in c(1,2,3,4,5)) cat(sprintf("     x = %2.0f : log((1-n)/n) - x^2 = %.10f\n", xx, lg(xx)-xx^2))

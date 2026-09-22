# Does the mass ceiling depend on trusting the state at arbitrarily high momentum?
#
# A referee's question the paper does not currently answer. I = (1/pi^2) int x^2 n(x) dx runs to
# infinity, and M_1 goes as I^(-2/5). If most of I came from the far ultraviolet, the bound would
# be hostage to physics nobody has. This is a cutoff-sensitivity check in the renormalisation-group
# spirit: integrate out above X and see what moves. It is NOT an RG calculation; there is no
# running coupling and nothing is resummed.
nmin <- function(x) (1-sqrt(1-exp(-x^2)))/2
Ifun <- function(hi) integrate(function(x) x^2*nmin(x), 0, hi, subdivisions=8000)$value/pi^2
I0   <- Ifun(Inf); M0 <- 491.6
cat(sprintf("  full integral I0 = %.10f, endpoint %.1f PeV\n\n", I0, M0))
cat("      cutoff X    I(<X)/I0     mass if truncated    shift from 491.6\n")
for (X in c(1,1.5,2,2.5,3,4,5,6)) {
  f <- Ifun(X)/I0
  cat(sprintf("   %10.1f  %10.7f  %17.2f  %15.3f PeV\n", X, f, M0*f^(-2/5), M0*f^(-2/5)-M0))
}
cat("\n  The integrand peaks where the band is wide and the measure is growing; locate it:\n")
g <- function(x) x^2*nmin(x)
xs <- seq(0.01,6,length.out=4000); pk <- xs[which.max(g(xs))]
cat(sprintf("     x^2 n(x) peaks at x = %.3f ; half the integral is below x = %.3f\n",
    pk, uniroot(function(v) Ifun(v)/I0-0.5, c(0.1,5))$root))
cat(sprintf("     99 per cent below x = %.3f ; 99.9 per cent below x = %.3f\n",
    uniroot(function(v) Ifun(v)/I0-0.99, c(0.5,8))$root,
    uniroot(function(v) Ifun(v)/I0-0.999, c(0.5,10))$root))
cat("\n=== which moment does the ultraviolet condition actually police? ===\n")
# m = 2 is the production integral; m = 3 is the energy. Compare the adopted Gaussian state
# against a p^-4 tail, on both moments, as a function of cutoff.
tail <- function(x) 0.5/(1+(x/0.62)^4)
mom  <- function(f,m,hi) integrate(function(x) x^m*f(x), 1e-6, hi, subdivisions=8000)$value
cat("      cutoff X   m=2 Gaussian   m=2 p^-4     m=3 Gaussian   m=3 p^-4\n")
for (X in c(3,10,30,100,300)) {
  cat(sprintf("   %10.0f  %12.6f %11.4f  %14.6f %10.4f\n",
      X, mom(nmin,2,X), mom(tail,2,X), mom(nmin,3,X), mom(tail,3,X)))
}
cat("\n  READING, stated with the moments separated because an earlier draft of this file ran\n")
cat("  them together. The PRODUCTION integral (m = 2) converges for BOTH states: x^2 * x^-4 is\n")
cat("  x^-2, which is integrable, so a p^-4 tail does not diverge there. What the ultraviolet\n")
cat("  condition polices is the ENERGY (m = 3): x^3 * x^-4 is x^-1, and the p^-4 column above\n")
cat("  grows without bound as the cutoff rises while the Gaussian column stops moving.\n")
cat("\n  So the two facts are about different moments and both hold. The number this paper quotes\n")
cat("  is an infrared quantity: half of I comes from below x = 1.03, 99.9 per cent from below\n")
cat("  x = 2.83, and truncating above x = 3 moves the endpoint by 0.076 PeV, twenty-six times\n")
cat("  smaller than its own quoted width. The ceiling therefore does not depend on the state at\n")
cat("  large momentum. The ultraviolet condition is not idle either, but it earns its keep on\n")
cat("  the energy rather than on the abundance.\n")

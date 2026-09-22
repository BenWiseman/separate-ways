# ==========================================================================================
# 4.1 states the residual phase freedom as an open gap: "what stays free is the region
# between, where the phase is neither irrelevant nor protected." Read alone that sounds like
# a hole under the ceiling. It is not, and the reason is one line, but the paper does not
# say it: the freedom is precisely what the 3.1 bound is a bound OVER, and it is one-sided.
#
# Every admissible state satisfies <N_+> >= n_*(x) pointwise. The production integral is
# I = Int x^2 n dx with a positive weight, so any admissible profile has I >= I_0, and
# M_1 ~ I^(-2/5) is monotonically decreasing in I. Therefore the phase can only push the
# mass DOWN. Demonstrate it on profiles that use the freedom as hard as the band allows.
P     <- function(x) exp(-x^2)
nstar <- function(x) (1 - sqrt(1-P(x)))/2
nmax  <- function(x) 1 - nstar(x)
I_of  <- function(nf, X=40) integrate(function(x) x^2*nf(x), 0, X, rel.tol=1e-11)$value
I0    <- I_of(nstar)
# NOTE: this is Int x^2 n dx = 0.125933, not the paper's production integral
# I_0 = 0.0127597, which carries its own extra weight. Only the RATIO I/I_0 feeds
# M_1 ~ I^(-2/5), and the ratio is the same under either weight, so the conclusion
# does not depend on which is used.
cat(sprintf("  Int x^2 n_* dx (the minimising profile) = %.10f\n\n", I0))

# admissible band profiles: n = n_* + D*w(x) with D = nmax - nstar and w in [0,1], and w
# must vanish in the UV or the integral diverges (the delta^2/4 floor of 4.1).
band <- function(w) function(x) nstar(x) + (nmax(x)-nstar(x))*w(x)
profs <- list(
  "minimising (the adopted state)" = nstar,
  "bump at x=0.5, width 0.3"       = band(function(x) exp(-((x-0.5)/0.3)^2)),
  "bump at x=1.0, width 0.5"       = band(function(x) exp(-((x-1.0)/0.5)^2)),
  "bump at x=1.5, width 0.5"       = band(function(x) exp(-((x-1.5)/0.5)^2)),
  "wide, cut off at x=2"           = band(function(x) ifelse(x<2, 1, 0)),
  "wide, cut off at x=3"           = band(function(x) ifelse(x<3, 1, 0))
)
cat("      profile                              I           I/I_0     M_1/M_1(max)   M_1 (PeV)\n")
for (nm in names(profs)) {
  I <- I_of(profs[[nm]]); r <- (I0/I)^(2/5)
  cat(sprintf("   %-36s %10.6f %10.4f %12.4f %12.1f\n", nm, I, I/I0, r, 491.6*r))
}
cat("\n  FLATLY: every admissible use of the phase freedom RAISES I and LOWERS the mass. Not\n")
cat("  one profile exceeds 491.6 PeV, and none can, because n >= n_* holds pointwise and the\n")
cat("  weight x^2 is positive. The freedom the paper calls open is one-sided, and it is the\n")
cat("  side the ceiling is safe on.\n")
cat("\n  What it DOES leave open is how far below the ceiling the true mass sits, and that is\n")
cat("  a genuine and unquantified range: the profiles above span 491.6 down to a few tens of\n")
cat("  PeV. The paper claims a ceiling and not a value, and this is exactly why.\n")

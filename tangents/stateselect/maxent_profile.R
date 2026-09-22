# Converging lines: an independent cross-field review proposes MaxEnt with moment constraints for
# the free profile, and section 6 already notes a Gaussian is the least-committed choice at fixed
# second moment. Test whether that is real or a coincidence of wording.
#
# The occupation is not a distribution, but x^2 n(x)/I IS one: the normalised distribution of
# produced particles over momentum, which is what the production integral integrates. Ask what
# distribution that is.

n  <- function(x) (1-sqrt(1-exp(-x^2)))/2
I  <- integrate(function(x) x^2*n(x), 0, Inf)$value
f  <- function(x) x^2*n(x)/I                      # normalised momentum distribution
m2 <- integrate(function(x) x^2*f(x), 0, Inf)$value
cat(sprintf("   normalisation  int f dx = %.10f\n", integrate(f,0,Inf)$value))
cat(sprintf("   second moment  <x^2>    = %.6f\n\n", m2))

# Maxwell (= MaxEnt on 3D momentum at fixed mean energy) with the SAME second moment
a  <- sqrt(m2/3)
mw <- function(x) sqrt(2/pi)*x^2*exp(-x^2/(2*a^2))/a^3
cat(sprintf("   Maxwell with the same <x^2>: a = %.6f, normalisation %.8f\n",
    a, integrate(mw,0,Inf)$value))

xs <- seq(0.001, 8, length.out=4000)
d  <- f(xs); m <- mw(xs)
cat(sprintf("   max |f - Maxwell|      = %.5f   (peak of f is %.5f)\n", max(abs(d-m)), max(d)))
cat(sprintf("   L1 distance            = %.5f\n", sum(abs(d-m))*(xs[2]-xs[1])))
kl <- sum(ifelse(d>0 & m>0, d*log(d/m), 0))*(xs[2]-xs[1])
cat(sprintf("   KL(f || Maxwell)       = %.5f nats\n", kl))

cat("\n=== what this does and does not show\n\n")
cat("  The produced-particle momentum distribution is CLOSE to Maxwellian but not equal to it.\n")
cat("  It cannot be equal: n(0) = 1/2 exactly, forced by the contact condition, whereas a\n")
cat("  Maxwell profile has n falling from the start. The deviation is concentrated at small x,\n")
cat("  which is exactly where the fold pins the occupation and thermodynamics does not.\n\n")
cat("  So the MaxEnt route does NOT close the remaining freedom by itself. Maximum entropy at\n")
cat("  fixed second moment returns the Maxwell form; the actual distribution departs from it in\n")
cat("  the infrared by construction, because Theta-invariance fixes n(0) = 1/2 there and no\n")
cat("  entropy principle produced that. The two agree where neither is doing any work.\n")
cat("  Recorded as a negative so the route is not re-run: the convergence between an outside\n")
cat("  reviewer's suggestion and our own note is a convergence of LANGUAGE, not of result.\n")

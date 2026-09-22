# ==========================================================================================
# THE SEPARABLE CEILING, DERIVED WITHOUT AN EXTREME-POINT ARGUMENT.
#
# 3.1 previously justified the separable competitor by "a linear functional is minimised at an
# extreme point, so separability reduces to product pure states". A pre-submission review
# objected, correctly: intersecting the separable set with the CONTACT CONSTRAINT changes the
# feasible set, and the extreme points of an intersection need not be the extreme points of
# either factor. The number survives; the argument did not.
#
# THE CONSTRAINED DERIVATION, which uses only two facts and no geometry of the separable set.
#   (i)  Separability bounds the even-block coherence:  |rho_{00,11}| <= n(1-n).
#   (ii) The contact condition <N_+> = <N_-> forces    Re rho_{00,11} = P(n - 1/2)/sqrt(P(1-P)).
# Saturating (i) with (ii) gives the separable floor as the root of
#       P(1/2 - n)/sqrt(P(1-P)) = n(1-n).
# ==========================================================================================
nsep <- function(P) {
  if (P < 1e-8) return(sqrt(P)/2)   # as P -> 0 the constraint degenerates; n_sep -> sqrt(P)/2
  g <- function(n) P*(0.5-n)/sqrt(P*(1-P)) - n*(1-n)
  if (g(1e-14)*g(0.5) > 0) return(sqrt(P)/2)
  uniroot(g, c(1e-14, 0.5), tol=1e-15)$root
}
cat("  The root, with both sides of the saturated constraint printed so the solve is visible:\n\n")
cat("        P        n_sep         n_*         P(1/2-n)/sqrt(P(1-P))    n(1-n)\n")
for (P in c(0.5, 0.2, 0.05)) {
  n <- nsep(P)
  cat(sprintf("   %7.2f %12.9f %11.9f %20.9f %12.9f\n", P, n, (1-sqrt(1-P))/2,
              P*(0.5-n)/sqrt(P*(1-P)), n*(1-n)))
}
cat("\n  At P = 0.2 the root is 0.190983006, which is the value the review reached\n")
cat("  independently from the same separability bound.\n\n")
cat("  Small-P branch, where the constraint degenerates and n_sep -> sqrt(P)/2:\n")
for (P in c(1e-2,1e-4,1e-6))
  cat(sprintf("    P = %.0e   n_sep = %.10f   sqrt(P)/2 = %.10f\n", P, nsep(P), sqrt(P)/2))
Px <- function(x) exp(-x^2)
I0  <- integrate(function(x) x^2*(1-sqrt(1-Px(x)))/2, 0, 40, rel.tol=1e-12)$value
Isep<- integrate(function(x) sapply(x, function(xx) xx^2*nsep(Px(xx))), 0, 40,
                 rel.tol=1e-10, subdivisions=4000)$value
cat(sprintf("\n  Int x^2 n_*    = %.12f\n  Int x^2 n_sep  = %.12f\n  ratio          = %.6f\n",
            I0, Isep, Isep/I0))
cat(sprintf("  mass ratio (I0/Isep)^(2/5) = %.6f  ->  %.4f PeV\n",
            (I0/Isep)^(2/5), 491.6*(I0/Isep)^(2/5)))
# 3.1 quotes these divided by pi^2, which is the production integral's own convention
# (I = (1/pi^2) Int x^2 n dx). Print both so a reader comparing the file with the text is
# not left converting. The 0.0569210 the paper states is Isep/pi^2 and appears nowhere above.
cat(sprintf("  in the paper's convention, I = (1/pi^2) Int x^2 n dx:\n"))
cat(sprintf("     I_0   = %.7f   (3.1 and 2.2 quote 0.0127597)\n", I0/pi^2))
cat(sprintf("     I_sep = %.7f   (3.1 quotes 0.0569210)\n", Isep/pi^2))
stopifnot(abs(I0/pi^2 - 0.0127597) < 5e-7, abs(Isep/pi^2 - 0.0569210) < 5e-7)
cat("\n  Same 4.461007 and same 270.295 PeV as before, with the extreme-point step removed.\n")
cat("\n  A NOTE ON WHAT NOT TO DO HERE. A direct numerical minimisation over separable\n")
cat("  contact-satisfying states with a penalty method and Nelder-Mead returned 0.3600 at\n")
cat("  P = 0.2, against the true 0.1910: it found a local optimum and reported it as the\n")
cat("  minimum. The analytic route above is used because the numerical one was not\n")
cat("  trustworthy, and that is worth recording rather than hiding.\n")

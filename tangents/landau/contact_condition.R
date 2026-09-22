# The cheap, robust half of the power-law question, separated out so it does not wait on the
# expensive transmission fit.
#
# CLAIM: the half-and-half condition Theta imposes at a contact does not depend on the
# geometry there. It needs only that H is real, the diagonal odd across the contact and the
# off-diagonal even. Both of the fold's contacts satisfy that: the bang, where the diagonal
# is gamma*eta, and the singular locus A.15 identifies, where a Kasner interior makes it a
# power law in proper time. Check the matrix identity over a range of exponents.

sx <- matrix(c(0,1,1,0),2,2)
cat("   n       max |H(-t) - sx H(t) sx|     diagonal odd?    condition holds?\n")
for (n in c(0.25,0.5,1,1.5,2,3,4)) {
  worst <- 0
  for (t in c(0.2,0.7,1.3,2.5)) for (p in c(0.3,0.9)) {
    D <- sign(t)*abs(t)^n; Dm <- sign(-t)*abs(-t)^n
    H <- matrix(c(D,p,p,-D),2,2); Hm <- matrix(c(Dm,p,p,-Dm),2,2)
    worst <- max(worst, max(abs(Hm - sx%*%H%*%sx)))
  }
  cat(sprintf("  %5.2f %26.1e %17s %20s\n", n, worst, "yes", if (worst<1e-14) "YES" else "no"))
}
cat("\n  So Theta: psi(t) -> sx psi*(-t) preserves the equation for every exponent, Theta^2 = 1,\n")
cat("  and ray-invariance forces |psi_1(0)| = |psi_2(0)| = 1/sqrt2 at the contact regardless of\n")
cat("  what the geometry is doing there. One condition, both contacts, no geometry input.\n")

cat("\n=== what the sweep exponent costs, if a contact is not radiation-linear\n\n")
cat("  Scaling: the crossing time is t_c ~ (p/A)^(1/n) and the energy there is p, so the\n")
cat("  adiabatic parameter goes as p^((n+1)/n) and P ~ exp(-c p^q) with q = (n+1)/n.\n")
cat("  Every q > 0 is super-polynomial, so Hadamard passes at EVERY moment for EVERY n.\n")
cat("  FLAT NEGATIVE: the ultraviolet condition does not select the exponent, so it does not\n")
cat("  pick radiation out. What it does do is fix the SHAPE once the exponent is known.\n\n")
nmin <- function(x,q) (1-sqrt(1-exp(-x^q)))/2
I0 <- (1/pi^2)*integrate(function(x) x^2*nmin(x,2), 0, Inf)$value
cat(sprintf("   reference: q = 2 gives I = %.7f, the adopted 0.0127597\n\n", I0))
cat("      n      q=(n+1)/n       I(q)      I/I(2)    M_1 factor    M_1 (PeV)   line (PeV)\n")
for (n in c(0.5,0.75,1,1.5,2,3,6)) {
  q <- (n+1)/n
  I <- (1/pi^2)*integrate(function(x) x^2*nmin(x,q), 0, Inf)$value
  f <- (I/I0)^(-2/5)
  cat(sprintf("  %5.2f %12.3f %12.7f %10.4f %13.4f %12.1f %12.1f\n",
      n, q, I, I/I0, f, 491.6*f, 245.8*f))
}
cat("\n  Reading: n = 1 is radiation, gives q = 2, and returns the adopted Gaussian and the\n")
cat("  adopted numbers. So the Gaussian occupation is not an imported profile at all: it is\n")
cat("  what a LINEAR crossing gives, and the crossing is linear because the bang is radiation\n")
cat("  dominated. A.5 already argues radiation is the saturating case at the bang on separate\n")
cat("  grounds; this says what the choice is worth in PeV.\n")

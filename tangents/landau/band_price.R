# What the remaining import is WORTH, in PeV.
#
# From theta_family.R the Theta-invariant family spans, at each p,
#     n in [n_min, n_max],  n_min = (1-sqrt(1-P))/2,  n_max = (1+sqrt(1-P))/2,
# and the run showed n_min + n_max = 1 to seven figures, which is unitarity:
# (|a1|-|a2|)^2/2 + (|a1|+|a2|)^2/2 = |a1|^2+|a2|^2 = 1.
# So the band is symmetric about 1/2. Two ends are already pinned:
#   deep IR: at p=0 the Hamiltonian is diagonal, no mixing, n = 1/2 for EVERY phase,
#            so the band has zero width there;
#   deep UV: Hadamard forces n -> n_min (mu_uv.R).
# Only the middle is free. Price the worst excursion.

nmin <- function(x) (1-sqrt(1-exp(-x^2)))/2
nmax <- function(x) (1+sqrt(1-exp(-x^2)))/2
Iof  <- function(f) (1/pi^2)*integrate(function(x) x^2*f(x), 0, Inf, subdivisions=4000)$value

Imin <- Iof(nmin)
cat(sprintf("   I at the least-occupied member : %.7f   (5.1 uses 0.0127597)\n", Imin))
cat(sprintf("   I at the most-occupied member  : divergent, since n_max -> 1 and int x^2 dx diverges\n"))
cat("   So the upper branch is not admissible at all: it has no finite number density.\n")
cat("   The excursion must therefore be confined to a band. Price it.\n\n")

cat("   Worst case: sit at n_max out to x_c, then drop to n_min beyond.\n\n")
cat("      x_c      I(x_c)     I/I_min    M_1 factor    M_1 (PeV)    nu line (PeV)\n")
for (xc in c(0.10,0.25,0.50,0.75,1.00,1.25,1.50,2.00,3.00)) {
  a <- (1/pi^2)*integrate(function(x) x^2*nmax(x), 0, xc)$value
  b <- (1/pi^2)*integrate(function(x) x^2*nmin(x), xc, Inf, subdivisions=4000)$value
  I <- a+b; fac <- (I/Imin)^(-2/5)
  cat(sprintf("   %6.2f %11.6f %11.4f %12.4f %12.1f %14.1f\n",
      xc, I, I/Imin, fac, 491.6*fac, 245.8*fac))
}

cat("\n   And the same table read as the question a referee will ask:\n")
cat("   how far can the mass fall before the KM3NeT comparison stops working?\n\n")
KM <- 220
for (xc in c(0.25,0.50,0.75,1.00)) {
  a <- (1/pi^2)*integrate(function(x) x^2*nmax(x), 0, xc)$value
  b <- (1/pi^2)*integrate(function(x) x^2*nmin(x), xc, Inf, subdivisions=4000)$value
  fac <- ((a+b)/Imin)^(-2/5)
  cat(sprintf("   x_c = %.2f -> line at %.1f PeV, KM3NeT median %d PeV, %s\n",
      xc, 245.8*fac, KM, if (245.8*fac >= KM) "still above the median" else "BELOW the median: comparison breaks"))
}

cat("\n=== the direction that matters, stated plainly\n\n")
cat("  Every excursion INCREASES I and so DECREASES the mass. There is no way to push the\n")
cat("  mass up. So 491.6 PeV is the endpoint of the family the FOLD defines, not merely of\n")
cat("  the family BFT parametrise, and the remaining import can only move the number one way.\n")
cat("  That is worth having: an upper bound that the symmetry supplies is a different object\n")
cat("  from a preferred value that a prescription supplies.\n")

cat("\n=== the inversion: KM3NeT already bounds the remaining freedom\n\n")
nminf <- function(x) (1-sqrt(1-exp(-x^2)))/2
nmaxf <- function(x) (1+sqrt(1-exp(-x^2)))/2
Im <- (1/pi^2)*integrate(function(x) x^2*nminf(x), 0, Inf, subdivisions=4000)$value
line <- function(xc) {
  a <- (1/pi^2)*integrate(function(x) x^2*nmaxf(x), 0, xc)$value
  b <- (1/pi^2)*integrate(function(x) x^2*nminf(x), xc, Inf, subdivisions=4000)$value
  245.8*((a+b)/Im)^(-2/5)
}
xstar <- uniroot(function(xc) line(xc)-220, c(0.05, 2), tol=1e-10)$root
cat(sprintf("   The endpoint falls to KM3NeT's 220 PeV median at x_c = %.4f.\n", xstar))
cat(sprintf("   Check: line(x_c) = %.2f PeV.\n", line(xstar)))
cat("\n   Read it forwards: the endpoint is an UPPER bound, and the reconstructed 220 PeV\n")
cat("   event sits inside it. Read it backwards: any state deviating maximally from the\n")
cat(sprintf("   least-occupied member out to x_c > %.2f pushes the endpoint BELOW the event,\n", xstar))
cat("   so the event would exceed the model's own bound. The measurement therefore limits\n")
cat("   how far the surviving phase freedom can run, without any appeal to minimum energy.\n")
cat("   The import is not merely smaller than it was. It is now bounded by an observation.\n\n")
cat("   Caveats kept in view: the 220 PeV figure is a reconstructed median assuming an\n")
cat("   incident E^-2 spectrum, the worst-case profile used here sits at the band edge and\n")
cat("   no real state need do that, and this is the two-level mode statement of A.1.\n")

# Ben has raised the quantum Zeno effect repeatedly: accumulated unobserved measurement across a
# universe causing a misalignment. The sharp version here: weak measurement decoheres the PHASE
# mu, and mu is exactly what picks the adopted state out of the Theta-invariant family. So ask
# what phase decoherence does to the occupation. Compute, do not hope.
#
# n(mu) = |a1 + e^{i mu} a2|^2/2 = (1 + 2|a1||a2| cos(mu - phi))/2, and unitarity gave
# |a1|^2+|a2|^2 = 1, so the band is symmetric about 1/2 with half-width |a1||a2|.
# The minimum is n_min = (1 - sqrt(1-P))/2 with P = exp(-x^2), so 2|a1||a2| = sqrt(1-P).

P  <- function(x) exp(-x^2)
nm <- function(x) (1-sqrt(1-P(x)))/2
hw <- function(x) sqrt(1-P(x))/2          # band half-width

cat("=== 1. the phase-averaged state\n\n")
cat("   Averaging n(mu) over mu uniformly kills the cosine, leaving <n> = 1/2 at EVERY x.\n")
cat("   Check numerically against the band:\n\n")
cat("      x      n_min        n_max      (n_min+n_max)/2    <n>_mu by quadrature\n")
for (x in c(0.2,0.6,1.0,2.0,4.0)) {
  a <- integrate(function(m) (1+2*hw(x)*cos(m))/2, 0, 2*pi)$value/(2*pi)
  cat(sprintf("  %5.1f %11.6f %12.6f %16.6f %20.6f\n", x, nm(x), 1-nm(x), 0.5, a))
}
cat("\n   So COMPLETE phase decoherence sends the occupation to 1/2 at every momentum. The\n")
cat("   number density int x^2 n dx and the energy int x^3 n dx both diverge. A fully\n")
cat("   dephased state is not admissible at all.\n")

cat("\n=== 2. partial decoherence: how precise must the phase be?\n\n")
cat("   Near the minimum, n(mu_* + d) = n_min + sqrt(1-P) (1-cos d)/2 ~ n_min + d^2/4 at large x.\n")
cat("   The Gaussian tail is n_min ~ exp(-x^2)/4, so a phase error d puts a FLOOR of d^2/4 on\n")
cat("   the occupation, and the tail stops being Gaussian where the floor overtakes it:\n")
cat("       exp(-x^2)/4 < d^2/4   i.e.   x > sqrt(2 ln(1/d)).\n\n")
cat("      phase error d      x_max = sqrt(2 ln(1/d))     floor on n      tail beyond x_max\n")
for (d in c(1e-2,1e-5,1e-10,1e-20,1e-40)) {
  xm <- sqrt(2*log(1/d))
  cat(sprintf("   %14.0e %24.3f %16.2e %22s\n", d, xm, d^2/4, "constant, not Gaussian"))
}
cat("\n   The requirement is exponential in x: coherence to exp(-x^2/2). Even a phase good to\n")
cat("   one part in 1e40 only buys a Hadamard tail out to x ~ 13.6. Beyond any such x the\n")
cat("   occupation has a constant floor and the energy integral diverges like x^4.\n")

cat("\n=== 3. is that fatal? No, and the reason is specific\n\n")
cat("   The state is pure at the bang by construction, and the field is free afterwards, so the\n")
cat("   occupations are fixed at production and later decoherence does not revisit them. The\n")
cat("   decoherence 4.2 computes is between the expanding and contracting BRANCHES at\n")
cat("   t_dec ~ 1.4e-32 s, which is a different object from a mode's relative phase.\n")
cat("   So the Zeno worry is real in principle and is answered by purity at the contact, not\n")
cat("   by the effect being small. If anything interacted with these modes at the bang, the\n")
cat("   model's ultraviolet behaviour would fail immediately rather than gradually.\n")

cat("\n=== 4. and it identifies what the leftover freedom IS\n\n")
cat("      x        band half-width     n_min          does the phase matter here?\n")
for (x in c(0.05,0.2,0.5,0.9682,2,4)) {
  w <- hw(x)
  tag <- if (w < 0.05) "no: band is collapsed" else if (nm(x) < 1e-3) "no: Hadamard pins it" else "YES"
  cat(sprintf("  %6.3f %16.6f %14.3e %32s\n", x, w, nm(x), tag))
}
cat("\n   The phase is irrelevant in the infrared because the band has closed there (n = 1/2\n")
cat("   whatever mu is), and it is pinned in the ultraviolet because Hadamard forbids the\n")
cat("   floor. It matters only in between, which is the same finite band section 6 already\n")
cat("   names as the surviving import. The leftover freedom is exactly the region where the\n")
cat("   phase is neither irrelevant nor protected.\n")

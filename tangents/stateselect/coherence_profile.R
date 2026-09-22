# The last imported assumption is 5.1's minimum-occupation prescription. Section 6 says what
# is imported is "the infrared profile". That is looser than it needs to be, and rewriting
# the adopted occupation in the pairing's own variables shrinks it.
#
# BFT's occupation, in 4.2's dimensionless momentum x:
#     n(x) = (1 - sqrt(1 - exp(-x^2)))/2
# Solve it backwards. With 1-2n = sqrt(1-exp(-x^2)):
#     (1-2n)^2 = 1 - exp(-x^2)   =>   exp(-x^2) = 1 - (1-2n)^2 = 4n(1-n).
# So the prescription is not a statement about a level at all. It is

cat("   4 n(1-n) = exp(-x^2)     exactly, for every x.\n\n")

x  <- seq(0.01, 6, length.out = 20000)
n  <- (1 - sqrt(1 - exp(-x^2)))/2
lhs <- 4*n*(1-n)
rhs <- exp(-x^2)
cat(sprintf("   max |4n(1-n) - exp(-x^2)| over x in [0.01,6] = %.3e\n", max(abs(lhs-rhs))))

cat("\n=== what 4n(1-n) is, for a mode paired across the sheets\n\n")
cat("  A Gaussian two-sheet state is a two-mode squeezed pair. For a fermionic pair mode\n")
cat("  with occupation n on each sheet, the pure-state pair is\n")
cat("     |psi> = sqrt(1-n)|0,0> + sqrt(n)|1,1>,\n")
cat("  whose off-diagonal coherence between the sheets is sqrt(n(1-n)), so the squared\n")
cat("  coherence normalised to its maximum is exactly 4n(1-n). Checking that the algebra\n")
cat("  and the reduced density matrix agree:\n\n")
C2 <- 4*n*(1-n)
# reduced state of one sheet is diag(1-n, n); its purity is (1-n)^2+n^2 = 1 - 2n(1-n)
pur <- (1-n)^2 + n^2
cat(sprintf("   max |4n(1-n) - 2(1 - purity)| = %.3e   (identity: 1-purity = 2n(1-n))\n",
            max(abs(C2 - 2*(1-pur)))))

cat("\n=== so the imported prescription, restated\n\n")
cat("  'Minimise the occupation' is, mode by mode and with no approximation,\n\n")
cat("     the coherence between the two sheets is GAUSSIAN in momentum,\n")
cat("     C(p)^2 = exp(-x^2),   x = sqrt(pi) p / sqrt(gamma),   gamma = M_1 a_1.\n\n")
cat("  Two things follow that the 'infrared profile' phrasing hides.\n")
cat("  First, the SCALE is not imported. x is p/sqrt(gamma) and gamma = M_1 a_1 is fixed by\n")
cat("  the mass the abundance already determines, so what is imported is a dimensionless\n")
cat("  shape, not a scale.\n")
cat("  Second, the shape is a Gaussian, which is one of the two profiles that survive the\n")
cat("  Hadamard moment test at every m. The imported piece therefore does not merely\n")
cat("  happen to be consistent with the ultraviolet condition; it is in the small class\n")
cat("  that condition leaves.\n")

cat("\n=== does the moment condition pick the Gaussian out of that class? Test it.\n\n")
cat("  Compare the Gaussian against an exponential coherence with the same second moment,\n")
cat("  both of which pass every moment. If the condition cannot separate them, say so.\n\n")
m2g <- sum(x^2 * exp(-x^2)) / sum(exp(-x^2))
lam <- sqrt(2/m2g)                       # exp(-lam x) matched on <x^2>? solve numerically
f <- function(l) sum(x^2*exp(-l*x))/sum(exp(-l*x)) - m2g
l  <- uniroot(f, c(0.05, 50))$root
cat(sprintf("   Gaussian <x^2> = %.5f ;  exponential rate matched at lambda = %.5f\n", m2g, l))
cat("\n      m      sum n k^m  (Gaussian)     sum n k^m  (exponential)    both finite?\n")
for (mm in c(0, 2, 4, 8, 16)) {
  g <- sum(x^2 * x^mm * exp(-x^2))
  e <- sum(x^2 * x^mm * exp(-l*x))
  cat(sprintf("   %5d %22.4e %26.4e %14s\n", mm, g, e, "yes"))
}
cat("\n  FLAT NEGATIVE: the moment condition does not separate them. Both converge at every\n")
cat("  m, so Hadamard narrows the imported piece to a class and does not pick the member.\n")
cat("  That is the same verdict 5.1 reaches for the state and it is reached here for the\n")
cat("  profile, one level down. The import is smaller than section 6 says, and it is still\n")
cat("  an import.\n")

cat("\n=== the route that is left open, stated so it can be worked\n\n")
cat("  A Gaussian is the maximum-entropy profile at fixed second moment, and the second\n")
cat("  moment here is fixed by gamma, which the mass fixes. So the remaining import would\n")
cat("  be DERIVED, not assumed, by any argument that makes the pair state maximum-entropy\n")
cat("  at fixed gamma. That is a variational statement about inter-sheet entanglement, in\n")
cat("  the same currency A.13 works in, and it is the first version of this question that\n")
cat("  is posed in the fold's own variables rather than BFT's.\n")

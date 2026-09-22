# ==========================================================================================
# The Zeno question, asked properly. If the two sheets measure each other through the fold
# during the crossing, the pair block decoheres. Decoherence in a two-level crossing drives
# the populations toward equipartition, n -> 1/2, rather than toward the coherent n_*.
#
# Two things to settle, neither by assertion:
#   (1) Does decoherence break the 491.6 PeV ceiling? Which direction does it move the mass?
#   (2) Is there a limit in which it is not merely disfavoured but inconsistent?
# ==========================================================================================
P     <- function(x) exp(-x^2)
nstar <- function(x) (1 - sqrt(1-P(x)))/2

cat("  (1) DIRECTION. The operator inequality Q >= n_* 1 is a statement about a matrix, so it\n")
cat("      holds for EVERY state on the block, pure or mixed. A decohered state is a state on\n")
cat("      the block. So decoherence cannot push the occupation below n_*; it can only raise\n")
cat("      it toward 1/2. Check what that does to the mass, which goes as M ~ I^(-2/5).\n\n")
# partial dephasing: interpolate mode by mode toward equipartition
Ilam <- function(lam, X) integrate(function(x) x^2*((1-lam)*nstar(x) + lam/2), 0, X,
                                   rel.tol=1e-12)$value
I0 <- Ilam(0, 40)
cat(sprintf("      coherent I_0 (cutoff x<40) = %.10f\n", I0))
cat("      lambda    I(lambda)        I/I_0      M/M_0 = (I_0/I)^(2/5)\n")
for (lam in c(0, 1e-6, 1e-4, 1e-2, 0.1)) {
  I <- Ilam(lam, 40)
  cat(sprintf("   %8.0e %14.4f %12.2f %18.6f\n", lam, I, I/I0, (I0/I)^(2/5)))
}
cat("\n      Every entry has M/M_0 <= 1. Decoherence at the bang LOWERS the mass, so the\n")
cat("      ceiling is safe against it: 491.6 PeV stays an upper bound a fortiori.\n")

cat("\n  (2) CONSISTENCY. Now the cutoff dependence, which is the whole story.\n\n")
cat("      cutoff X     I(lambda=0)      I(lambda=1e-6)    I(lambda=1e-4)\n")
for (X in c(10, 20, 40, 80, 160)) {
  cat(sprintf("   %10.0f %16.8f %18.4f %18.4f\n", X, Ilam(0,X), Ilam(1e-6,X), Ilam(1e-4,X)))
}
cat("\n      The coherent integral converges. ANY constant nonzero lambda adds lambda*X^3/6,\n")
cat("      which diverges with the cutoff no matter how small lambda is. So a uniform\n")
cat("      decoherence rate across modes does not give a small correction to the dark-matter\n")
cat("      density; it gives an infinite one.\n")

cat("\n  So the requirement is on the SHAPE of lambda(x), not its size. For convergence,\n")
cat("  x^2 lambda(x) must be integrable, i.e. lambda must fall faster than x^-3. Check how\n")
cat("  much room the coherent solution leaves:\n\n")
cat("        x      n_*(x)        x^-3        n_*(x) / x^-3\n")
for (xx in c(2,3,4,5,6)) {
  cat(sprintf("   %6.1f %12.4e %12.4e %16.4e\n", xx, nstar(xx), xx^-3, nstar(xx)*xx^3))
}
cat("\n  FLATLY, and this is the part worth keeping: fold-invariance is not merely one way to\n")
cat("  make the production integral finite. Any decoherence at the bang that does not switch\n")
cat("  itself off in the ultraviolet faster than x^-3 makes the dark-matter density diverge.\n")
cat("  The Hadamard/ultraviolet criterion of 2.2 IS that requirement, arrived at from the\n")
cat("  other side. Coherence at the bang is what makes the answer finite, and the charge-odd\n")
cat("  decomposition says the same thing: the divergent piece is the vacuum half, which is\n")
cat("  charge-even, and the finite physical content is the charge-odd remainder.\n")

cat("\n  ==========================================================================\n")
cat("  WHICH DEPHASING PROFILES SURVIVE THE FINITENESS TEST?\n")
cat("  ==========================================================================\n")
cat("  The condition is x^2 lambda(x) integrable. Test candidate profiles rather than\n")
cat("  argue about them. C(x) = exp(-x^2/2) is the mode's own concurrence.\n\n")
Cx <- function(x) exp(-x^2/2)
profs <- list(
  "uniform, lambda = const"        = function(x) rep(1e-6, length(x)),
  "power law, lambda ~ x^-1"       = function(x) 1e-6*x^-1,
  "power law, lambda ~ x^-3"       = function(x) 1e-6*x^-3,
  "power law, lambda ~ x^-4"       = function(x) 1e-6*x^-4,
  "tracks coherence, lambda ~ C"   = function(x) 1e-6*Cx(x),
  "tracks coherence, lambda ~ C^2" = function(x) 1e-6*Cx(x)^2,
  "UV-growing, lambda ~ x^+1"      = function(x) 1e-6*x
)
cat(sprintf("   %-34s %14s %14s %14s   verdict\n", "profile", "X=20", "X=80", "X=320"))
for (nm in names(profs)) {
  f <- profs[[nm]]
  v <- sapply(c(20,80,320), function(X)
        integrate(function(x) x^2*f(x), 1e-3, X, rel.tol=1e-10, subdivisions=6000)$value)
  # x^-3 is the marginal case: Int x^2 x^-3 dx = Int dx/x = log x, so it diverges but
  # only logarithmically, and its column barely moves. Label it so nobody reads the slow
  # numbers as a contradiction of the verdict.
  ratio <- v[3]/v[1]
  verdict <- if (ratio < 1.05) "converges" else if (ratio < 2) "DIVERGES (log, marginal)" else "DIVERGES"
  cat(sprintf("   %-34s %14.3e %14.3e %14.3e   %s\n", nm, v[1], v[2], v[3], verdict))
}
cat("\n  So the dividing line is sharp and it is not about strength. A dephasing that tracks\n")
cat("  the coherence it destroys passes with room to spare, because a Gaussian beats every\n")
cat("  power. A dephasing with any power-law tail shallower than x^-3, and a fortiori one\n")
cat("  that grows in the ultraviolet, makes the dark-matter density infinite however weak\n")
cat("  it is made.\n")
cat("\n  WHAT THIS DOES NOT SETTLE, and it is the whole of the Zeno question: whether graviton\n")
cat("  exchange between the sheets has a rate that tracks the coherence or one that grows\n")
cat("  with momentum. Gravity couples more strongly at higher energy, which points the wrong\n")
cat("  way, but a bare coupling is not a renormalised decoherence rate and the Hadamard\n")
cat("  condition is precisely the statement that the ultraviolet structure is the flat-space\n")
cat("  one. Deciding it needs the interacting theory this programme does not yet have. The\n")
cat("  test is stated here so that whoever builds that theory has a number to hit: the\n")
cat("  inter-sheet dephasing must fall faster than the cube of the momentum.\n")

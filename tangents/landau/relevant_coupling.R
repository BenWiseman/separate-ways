# Two separate results in this paper say the same thing in different words:
#   - 4.1: a phase error delta puts a FLOOR delta^2/4 on the occupation, and the production
#     integral then diverges for any delta != 0.
#   - zeno_at_the_bang.R: any uniform inter-sheet dephasing lambda adds lambda X^3/6 and
#     diverges however small lambda is.
# Both are perturbations that cannot be small. That is the language of a RELEVANT coupling,
# and the renormalisation group has a name and a number for it. Extract the number.
P     <- function(x) exp(-x^2)
nstar <- function(x) (1 - sqrt(1-P(x)))/2
I_of  <- function(nf, X) integrate(function(x) x^2*nf(x), 0, X, rel.tol=1e-12)$value

cat("  How each contribution scales with the cutoff X. A term going as X^d has scaling\n")
cat("  dimension d: d > 0 is relevant, d = 0 marginal, d < 0 irrelevant.\n\n")
Xs <- c(10, 20, 40, 80)
cat("        cutoff X    coherent n_*      phase floor d^2/4   uniform dephasing\n")
for (X in Xs) {
  a <- I_of(nstar, X)
  b <- I_of(function(x) rep((1e-3)^2/4, length(x)), X)
  c_ <- I_of(function(x) rep(1e-6, length(x)), X)
  cat(sprintf("   %10.0f %16.10f %20.4f %18.4f\n", X, a, b, c_))
}
cat("\n  Fit the exponent from the last doubling:\n")
for (nm in c("coherent n_*","phase floor","dephasing")) {
  f <- switch(nm, "coherent n_*"=nstar,
                  "phase floor"=function(x) rep((1e-3)^2/4, length(x)),
                  "dephasing"=function(x) rep(1e-6, length(x)))
  r <- I_of(f, 80)/I_of(f, 40)
  cat(sprintf("   %-16s I(80)/I(40) = %10.4f  ->  d = %.3f\n", nm, r, log2(r)))
}
cat("\n  So the coherent occupation is irrelevant in the technical sense: its integral stops\n")
cat("  growing, d = 0 with a convergent coefficient. Both contaminations are RELEVANT with\n")
cat("  dimension 3, which in this variable is the volume factor x^2 dx integrated against a\n")
cat("  constant. Dimension 3 is strongly relevant: it cannot be made small by taking the\n")
cat("  coupling small, only by setting it to zero.\n")
cat("\n  WHAT THAT MEANS, and it is the part worth keeping. A relevant coupling tuned to zero\n")
cat("  is normally a fine-tuning problem. Here it is not, because it is tuned by a SYMMETRY:\n")
cat("  Theta-invariance sets the phase error to zero exactly, and does so for every mode at\n")
cat("  once. The fold is therefore doing the job a symmetry exists to do, protecting a\n")
cat("  strongly relevant operator whose presence would make the dark-matter density\n")
cat("  infinite rather than merely wrong.\n")
cat("\n  WHAT IT DOES NOT MEAN. This is a statement about the free two-level reduction and its\n")
cat("  cutoff behaviour, not a renormalisation-group analysis of an interacting theory. No\n")
cat("  beta function has been computed and no operator mixing has been considered. It names\n")
cat("  the structure that two existing results share; it does not add a new one.\n")

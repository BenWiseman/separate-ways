# Ben has asked repeatedly about MISALIGNMENT between the sheets - the quantum Zeno
# tangent is a misalignment question, and so is "do the sheets carry their own field".
# The algebraic reciprocal law gives the right object to measure it with.
#
# If sigma is the fold's involutive automorphism, the state's failure to be its own fold
# image is the relative entropy S(omega || omega o sigma). Compute it and see what the
# involution forces.
#
# In finite dimensions with density matrix rho and sigma(rho) = S rho S, S^2 = 1:
#   Srel(rho || S rho S) = tr[rho log rho] - tr[rho log(S rho S)]
#                        = tr[(rho - S rho S) log rho]

set.seed(71); n <- 10
mklog <- function(M) { e <- eigen(M, symmetric=TRUE); e$vectors %*% diag(log(pmax(e$values,1e-300))) %*% t(e$vectors) }
Srel  <- function(a,b) sum(diag(a %*% (mklog(a) - mklog(b))))
mkrho <- function() { M <- matrix(rnorm(n*n),n,n); M <- M %*% t(M) + 0.2*diag(n); M/sum(diag(M)) }
mkinv <- function() { Q <- qr.Q(qr(matrix(rnorm(n*n),n,n))); k <- sample(1:(n-1),1)
                      Q %*% diag(c(rep(1,k), rep(-1,n-k))) %*% t(Q) }

cat("=== 1. relative entropy is famously ASYMMETRIC. Is it here?\n\n")
cat("        trial   S(rho || s.rho)    S(s.rho || rho)      difference\n")
for (i in 1:6) {
  rho <- mkrho(); S <- mkinv(); sr <- S %*% rho %*% S
  a <- Srel(rho, sr); b <- Srel(sr, rho)
  cat(sprintf("   %8d %17.10f %18.10f %15.2e\n", i, a, b, abs(a-b)))
}
e <- max(replicate(400, { rho <- mkrho(); S <- mkinv(); sr <- S %*% rho %*% S
      abs(Srel(rho,sr) - Srel(sr,rho)) }))
cat(sprintf("\n   max asymmetry over 400 random (state, involution) pairs = %.2e\n", e))
cat("\n  Symmetric, exactly. The reason is two lines: relative entropy is invariant under\n")
cat("  a common automorphism, so S(s.rho || rho) = S(s.rho || s.(s.rho)) = S(rho || s.rho)\n")
cat("  once sigma^2 = id. Asymmetry is the normal state of affairs for relative entropy,\n")
cat("  so this is a property the fold confers rather than one it inherits.\n")

cat("\n=== 2. does it vanish exactly when the state is fold-symmetric?\n\n")
rho0 <- mkrho(); S <- mkinv()
sym <- (rho0 + S %*% rho0 %*% S)/2                       # fold-symmetric by construction
cat(sprintf("   fold-symmetric state:  |s.rho - rho| = %.2e,  S = %.3e\n",
    max(abs(S %*% sym %*% S - sym)), Srel(sym, S %*% sym %*% S)))
cat("\n        mixing t      |s.rho - rho|        S(rho || s.rho)      S / |.|^2\n")
for (t in c(0, 0.01, 0.05, 0.2, 0.6, 1)) {
  r <- (1-t)*sym + t*rho0; sr <- S %*% r %*% S
  d <- max(abs(sr - r)); v <- Srel(r, sr)
  cat(sprintf("   %11.2f %18.3e %20.3e %16s\n", t, d, v, if (d>1e-9) sprintf("%.3f", v/d^2) else "-"))
}
cat("\n  Vanishes only at exact fold symmetry and grows QUADRATICALLY in the misalignment,\n")
cat("  with the ratio settling. So a small failure of the fold costs second order, not\n")
cat("  first: the fold-symmetric states are a stationary point of this divergence, not\n")
cat("  merely its zero set.\n")

cat("\n=== 3. what that buys, and what it does not\n\n")
cat("  BUYS: a well-defined, symmetric, second-order measure of how far the universe\n")
cat("  falls short of being its own CPT image, which is the quantity Ben's misalignment\n")
cat("  questions are about. Second order matters: it means accumulating many small\n")
cat("  independent misalignments adds their SQUARES, so a universe-wide sum over N\n")
cat("  uncorrelated events grows like N rather than N^2, and per-event effects of size\n")
cat("  eps give a total of order N eps^2.\n\n")
cat("        N events        eps per event       total misalignment ~ N eps^2\n")
for (N in c(1e10, 1e20, 1e40)) for (eps in c(1e-20, 1e-30)) 
  cat(sprintf("   %14.0e %18.0e %30.2e\n", N, eps, N*eps^2))
cat("\n  DOES NOT BUY: any of those numbers. The per-event epsilon is not computed here\n")
cat("  and nothing in the fold supplies it, so the table is a scaling and not a result.\n")
cat("  Writing it as a prediction would be inventing the input.\n")

cat("\n=== 4. the route\n\n")
cat("  The quantity to compute is the per-event epsilon for a definite process, and the\n")
cat("  natural first one is the decoherence of section 4.1, where a branch that decoheres\n")
cat("  is exactly a state that stops being its own fold image. That is a calculation the\n")
cat("  paper already has the machinery for, and it would turn the scaling above into a\n")
cat("  number. Not attempted here.\n")

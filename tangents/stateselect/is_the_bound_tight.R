# A reviewer reading 3.1 as an SDP relaxation (NPA / Lasserre / DPS) will say the literature
# supplies a CONVERGENT HIERARCHY of stronger bounds. That is true of a relaxation. It is not
# true of a bound that is ATTAINED: no outer approximation improves on an infimum that some
# admissible state actually reaches. So the question is not "which level is this" but "is the
# minimiser admissible". Settle it by exhibiting the minimiser.
I2 <- diag(2); Z <- diag(c(1,-1)); f <- matrix(c(0,1,0,0),2,2,byrow=TRUE)
a <- kronecker(f, I2); b <- kronecker(Z, f); ad <- t(a); bd <- t(b)
Np <- (ad%*%a + bd%*%b)/2                      # per pair member, out region

cat("       P     min eig(Q)        n_*         <N_+>        <N_->     |<N_+>-<N_->|\n")
for (P in c(0.9,0.5,0.2,0.05,0.01)) {
  s <- sqrt(P); cc <- sqrt(1-P)
  am <- cc*a + s*bd; bm <- cc*b - s*ad
  Nm <- (t(am)%*%am + t(bm)%*%bm)/2            # per pair member, in region
  Q  <- (Np + Nm)/2
  e  <- eigen(Q, symmetric=TRUE)
  k  <- which.min(e$values); v <- e$vectors[,k]
  v  <- v/sqrt(sum(v^2))
  np <- as.numeric(t(v)%*%Np%*%v); nm <- as.numeric(t(v)%*%Nm%*%v)
  cat(sprintf("   %6.2f %12.9f %12.9f %12.9f %12.9f %14.2e\n",
              P, min(e$values), (1-sqrt(1-P))/2, np, nm, abs(np-nm)))
}
cat("\n  The minimising eigenvector satisfies the contact condition <N_+> = <N_-> to machine\n")
cat("  precision at every P, and its Q-expectation IS n_*. So the infimum is attained by an\n")
cat("  admissible state, the inequality Q >= n_* 1 is tight, and no relaxation hierarchy can\n")
cat("  return a larger infimum. Levels above the first would return the same number.\n")
cat("\n  What the SDP literature does give is a NAME and a proof style: this is a positivity\n")
cat("  certificate over a constrained state set, and the separable competitor of 3.1 is the\n")
cat("  DPS setting exactly, since the extreme points of the separable set are product pure\n")
cat("  states, which is why minimising a LINEAR functional there reduces to them. Useful for\n")
cat("  legibility to a quantum-information referee. Not a route to a stronger number.\n")

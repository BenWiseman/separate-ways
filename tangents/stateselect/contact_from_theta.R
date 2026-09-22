# ==========================================================================================
# An editorial pass objected that the contact condition <N_+> = <N_-> is "assumed, not
# derived", and that the sharpness of the ceiling rests on it. It follows from what the fold
# already is, in two lines, and the argument is analytic rather than numerical.
#
# THE ARGUMENT. Theta is antiunitary with Theta^2 = 1, and the fold exchanges the two
# asymptotic regions, so Theta N_- Theta^-1 = N_+. For a Theta-invariant state,
# Theta rho Theta^-1 = rho, and any hermitian A,
#
#     <Theta A Theta^-1> = Tr(rho Theta A Theta^-1)
#                        = Tr(Theta rho A Theta^-1)     (Theta rho Theta^-1 = rho)
#                        = Tr(rho A)^*                  (antilinearity conjugates the trace)
#                        = <A>                          (A hermitian, so <A> is real)
#
# With A = N_-, that is <N_+> = <N_->. The contact condition is a CONSEQUENCE of
# Theta-invariance, not an extra input.
#
# A first version of this file tried to "verify" that numerically by building an orthogonal
# map S with S N_- S^T = N_+ from eigenvectors and then checking S^2 = 1 and S N_- S^T = N_+.
# Both checks were worthless: the second is true by construction, and the first swung between
# 1e-16 and 2.0 because eigen() returns eigenvectors with arbitrary sign. Deleted. What
# follows is the one numerical statement that is not circular.
I2 <- diag(2); Z <- diag(c(1,-1)); f <- matrix(c(0,1,0,0),2,2,byrow=TRUE)
a <- kronecker(f,I2); b <- kronecker(Z,f); ad <- t(a); bd <- t(b)
Np <- (ad%*%a + bd%*%b)/2

cat("  The minimising state is obtained WITHOUT imposing the contact condition: it is the\n")
cat("  lowest eigenvector of Q. If the condition is a consequence rather than a selection,\n")
cat("  that state should satisfy it anyway. It does.\n\n")
cat("      P       <N_+>          <N_->        |difference|\n")
for (P in c(0.9,0.5,0.2,0.05,0.01,0.001)) {
  s <- sqrt(P); cc <- sqrt(1-P)
  am <- cc*a + s*bd; bm <- cc*b - s*ad
  Nm <- (t(am)%*%am + t(bm)%*%bm)/2
  Q  <- (Np + Nm)/2
  e  <- eigen(Q, symmetric=TRUE); v <- e$vectors[, which.min(e$values)]
  np <- as.numeric(t(v)%*%Np%*%v); nm <- as.numeric(t(v)%*%Nm%*%v)
  cat(sprintf("   %6.3f %13.9f %14.9f %16.2e\n", P, np, nm, abs(np-nm)))
}
cat("\n  WHAT IS STILL AN INPUT, and the paper should not claim otherwise: that the fold\n")
cat("  EXCHANGES the two asymptotic regions, Theta N_- Theta^-1 = N_+. That is what the fold\n")
cat("  is taken to be in 2.1, where the involution relates the sheets and the sheets are the\n")
cat("  two sides of the crossing. Given that, the contact condition follows. The derivation\n")
cat("  moves the assumption one step back, to the fold's own definition, which is where the\n")
cat("  paper wants its assumptions to live.\n")
cat("\n  So the tolerance analysis in 3.1 is not a patch on a guess. It prices what happens\n")
cat("  if Theta-invariance is only APPROXIMATE, which is a weaker and different question\n")
cat("  from whether the condition holds when Theta-invariance is exact.\n")

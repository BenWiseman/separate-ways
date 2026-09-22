# Check my own claim from the previous commit. I said the crossover shows "the two sheets are
# maximally entangled at long wavelength". That attribution needs testing, because at the BANG
# the two sheets are eta<0 and eta>0, which is a TIME direction: for a free field the data on
# one determines the data on the other, so they are one set of degrees of freedom at two times,
# not two subsystems you can trace between. Entanglement needs two commuting algebras.
#
# Build the state explicitly and see what the two entangled modes actually are.

# Fermionic pair production: beta mixes positive and negative frequency at fixed comoving p,
# producing a particle of momentum p and an antiparticle of momentum -p. The out-region state is
#    |out> = prod_p [ alpha_p |0_p 0_{-p}> + beta_p |1_p 1_{-p}> ],   |alpha|^2+|beta|^2 = 1.
# Trace over the partner and read the entropy.
nf   <- function(x) (1-sqrt(1-exp(-x^2)))/2
Sfun <- function(n) ifelse(n<=0|n>=1, 0, -n*log(n)-(1-n)*log(1-n))

cat("   x      n=|beta|^2   rho_1 = diag(1-n, n)   S(rho_1)    S(joint)   mutual info\n")
for (x in c(0.001, 0.5, 0.9682, 1.5, 3)) {
  n <- nf(x); a <- sqrt(1-n); b <- sqrt(n)
  psi <- c(a, 0, 0, b)                       # |00>, |10>, |01>, |11>
  rho <- outer(psi, Conj(psi))
  # partial trace over mode 2
  r1 <- matrix(c(rho[1,1]+rho[2,2], rho[1,3]+rho[2,4],
                 rho[3,1]+rho[4,2], rho[3,3]+rho[4,4]), 2, 2, byrow=TRUE)
  ev <- Re(eigen(r1)$values); ev <- ev[ev>1e-15]
  S1 <- -sum(ev*log(ev)); Sj <- 0
  cat(sprintf("  %6.3f %11.6f %22s %10.6f %10.6f %12.6f\n",
      x, n, sprintf("(%.4f, %.4f)", 1-n, n), S1, Sj, 2*S1))
  stopifnot(abs(S1 - Sfun(n)) < 1e-9)
}
cat("\n  S(rho_1) reproduces -n ln n - (1-n) ln(1-n) exactly at every x, and the joint state is\n")
cat("  pure, so the mutual information is 2S. The arithmetic of the crossover is unchanged.\n")

cat("\n=== but WHICH two modes are these?\n\n")
cat("  They are the two members of a produced PAIR: a particle of comoving momentum p and an\n")
cat("  antiparticle of momentum -p, both in the out region. They are two commuting subalgebras\n")
cat("  and tracing one out is legitimate.\n")
cat("  They are NOT the two sheets. At the bang the sheets are eta<0 and eta>0, related by\n")
cat("  evolution of one free field, so they are not independent subsystems and there is no\n")
cat("  entanglement entropy between them to compute. My previous commit attributed this\n")
cat("  entropy to the sheets. That was wrong and is corrected here.\n")

cat("\n=== what DOES survive, and it is most of it\n\n")
cat("  (1) The crossover itself: pair entanglement is maximal in the infrared, where n = 1/2\n")
cat("      exactly and S = ln 2, and collapses in the ultraviolet where n is Gaussian.\n")
cat(sprintf("      Half of the maximum is gone by x = %.4f, unchanged.\n",
      uniroot(function(x) Sfun(nf(x))-log(2)/2, c(0.1,4))$root))
cat("  (2) The scale coincidence: e^{-x^2} reaches 1/e at x = 1, within 3% of that point.\n")
cat("  (3) The SHEET statement that does hold is a symmetry, not an entanglement: Theta forces\n")
cat("      the occupation measured from the in-region to equal the one measured from the out-\n")
cat("      region, which is the n = m condition, and that is what picks the half-angle state.\n")
cat("  (4) At a HORIZON the two-subsystem reading is legitimate, because there J relates a\n")
cat("      wedge algebra to its commutant. The illegitimate step is importing that reading to\n")
cat("      the bang, where the two regions are not commuting algebras.\n")

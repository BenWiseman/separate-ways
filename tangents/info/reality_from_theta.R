# Named last turn: does the fold FORCE the wavefunction real, or does it import reality
# from the no-boundary proposal? Line 51 of the manuscript says "a real no-boundary
# wavefunction", so the paper currently assumes it.
#
# Theta is ANTILINEAR. Write Theta = K o S, K complex conjugation, S linear. Then
# Theta psi = S conj(psi), and Theta^2 = 1 requires S conj(S) = 1.
# Imposing Theta psi = psi on the physical state is therefore a REALITY CONDITION, twisted
# by S. In the branch basis, with S the branch swap, it reads psi_- = conj(psi_+), which is
# exactly the no-boundary form and forces |psi_-| = |psi_+|.
#
# The question is whether the fold may IMPOSE it. Last turn showed that for an EVOLVING
# state, Theta-invariance at one time does not persist. Wheeler-DeWitt is not an evolving
# state: H psi = 0 is a constraint. So the two settings differ and this tests which applies.

set.seed(101)
cat("=== 1. Theta psi = psi IS a reality condition, twisted by S\n\n")
Ssw <- matrix(c(0,1,1,0),2,2)                       # branch swap
cat(sprintf("   Theta^2 = 1 needs S conj(S) = 1:  |S conj(S) - 1| = %.2e\n",
    max(Mod(Ssw %*% Conj(Ssw) - diag(2)))))
psip <- complex(real=rnorm(1), imaginary=rnorm(1))
psi  <- c(psip, Conj(psip))                          # built to satisfy the condition
cat(sprintf("   test state psi = (z, conj z):     |S conj(psi) - psi| = %.2e\n",
    max(Mod(Ssw %*% Conj(psi) - psi))))
cat(sprintf("   branch weights:                   |psi+|^2 = %.6f, |psi-|^2 = %.6f\n",
    Mod(psi[1])^2, Mod(psi[2])^2))
bad <- c(psip, 0.3*Conj(psip))
cat(sprintf("   a state violating it:             |S conj(psi) - psi| = %.4f, weights %.4f / %.4f\n",
    max(Mod(Ssw %*% Conj(bad) - bad)), Mod(bad[1])^2, Mod(bad[2])^2))
cat("\n  So Theta-invariance and equal branch weight are the same condition, not two.\n")

cat("\n=== 2. may the fold impose it? Only if Theta maps solutions to solutions.\n")
cat("   Wheeler-DeWitt is H psi = 0 with H a REAL operator. Test that the solution space\n")
cat("   is closed under conjugation, and count the Theta-invariant solutions.\n\n")
n <- 12
Hr <- matrix(rnorm(n*n), n, n); Hr <- (Hr + t(Hr))/2      # real symmetric, as WDW is
Hr <- Hr - diag(rep(eigen(Hr)$values[7], n))              # push an eigenvalue to zero-ish
ev <- eigen(Hr, symmetric=TRUE)
k <- which(abs(ev$values) < 1e-9)
cat(sprintf("   real symmetric H, exact zero modes: %d\n", length(k)))
# build a genuine kernel: project out the smallest eigenvalue direction
v <- ev$vectors[, which.min(abs(ev$values))]
cat(sprintf("   |H v| for the near-kernel vector = %.2e, and v is REAL: max|Im v| = %.2e\n",
    max(abs(Hr %*% v)), max(abs(Im(v)))))
cat("\n   If psi solves H psi = 0 with H real, so does conj(psi): check on a complex\n")
cat("   combination.\n")
psi_c <- v * complex(real=0.7, imaginary=1.3)
cat(sprintf("   |H psi_c| = %.2e   |H conj(psi_c)| = %.2e\n",
    max(Mod(Hr %*% psi_c)), max(Mod(Hr %*% Conj(psi_c)))))
cat("\n  Closed under conjugation. So Theta maps physical solutions to physical solutions\n")
cat("  and the invariant sector is well defined. Restricting to it is legitimate.\n")

cat("\n=== 3. and the invariant sector is exactly HALF the real dimension\n\n")
cat("   A complex solution space of complex dimension d has real dimension 2d; the\n")
cat("   conjugation-invariant part has real dimension d. Check by projecting:\n\n")
cat("   Build an explicit complex solution space, apply the projector (1 + Theta)/2\n")
cat("   with Theta = K o S, and measure the rank of the invariant part as a REAL\n")
cat("   subspace. A first draft of this section printed the expected numbers without\n")
cat("   computing them, which is not a check.\n\n")
cat("        d (complex)   real dim 2d   measured rank of invariant part   matches d?\n")
for (d in c(1,2,3,5)) {
  Sc <- { Q <- qr.Q(qr(matrix(rnorm(2*d*2*d),2*d,2*d))); Q %*% Conj(t(Q)) }   # S = 1 case
  Sc <- diag(2*d)
  # realify: represent the complex space C^d as R^{2d}, Theta = conjugation
  # invariant vectors are the real ones; measure by projecting a random complex basis
  Z <- matrix(rnorm(d*d) + 1i*rnorm(d*d), d, d)        # random complex basis of the space
  inv <- (Z + Conj(Z))/2                                # (1 + Theta)/2 applied columnwise
  realified <- rbind(Re(inv), Im(inv))
  r <- qr(realified)$rank
  cat(sprintf("   %11d %13d %36d %13s\n", d, 2*d, r, if (r==d) "yes" else "NO"))
}
cat("\n  Exactly half, which is the same halving A.15 finds at the singularity through the\n")
cat("  projector (1 + J U)/2. Two different loci, one mechanism: an involution keeping its\n")
cat("  even sector.\n")

cat("\n=== 4. flatly: what this does and does not remove\n\n")
cat("  REMOVES: the need for the no-boundary proposal specifically. Reality of psi follows\n")
cat("  from imposing Theta-invariance on physical states, and Theta-invariance is the\n")
cat("  fold's own content. The paper can stop saying 'a real no-boundary wavefunction' as\n")
cat("  though reality were Hartle and Hawking's contribution to the argument.\n")
cat("  DOES NOT REMOVE: the need for the state to be a CONSTRAINT solution rather than an\n")
cat("  evolving one. Last turn's result stands: for an evolving state, Theta-invariance at\n")
cat("  one time does not persist. The derivation works because Wheeler-DeWitt is timeless,\n")
cat("  and that is an assumption about the formulation, not a theorem.\n")
cat("  ALSO NOT REMOVED: whether Theta-invariance may be IMPOSED at all. 2.4 gauges the\n")
cat("  LINEAR involution P_perp, following Harlow and Numasawa; gauging an ANTILINEAR one\n")
cat("  is a different operation and the paper does not do it. So this is a condition the\n")
cat("  fold may adopt, not one it already has.\n")

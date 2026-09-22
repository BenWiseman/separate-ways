# Named last turn: redo the branch imbalance with 4.1's ACTUAL setting rather than a
# random CPT-covariant matrix. Before computing anything new, reproduce 4.1's own two
# numbers, since one of them (the gauged threshold) was wrong in an earlier draft and an
# independent check is worth having.
#
# 4.1:  log|<E-|E+>| = -(N_f/4) sum_{n>=2} d_n log[1 + A^4(A^2-1)/(n^2(n^2-1)^2)],  A = aH
#   ungauged   d_n = n^2
#   gauged     d_n = sum over EVEN l of (2l+1), l = 0..n-1   (P_perp acts as (-1)^l on the
#              l multiplet, so it SPLITS each level rather than signing it)

dn_ungauged <- function(n) n^2
dn_gauged   <- function(n) { l <- 0:(n-1); sum((2*l+1)[l %% 2 == 0]) }

logov <- function(A, dfun, Nf=1, nmax=4000) {
  n <- 2:nmax
  d <- sapply(n, dfun)
  -(Nf/4)*sum(d * log(1 + A^4*(A^2-1)/(n^2*(n^2-1)^2))) }

cat("=== 0. the gauged dimension count, against the paper's stated n=2 case\n\n")
cat("        n     n^2 (ungauged)    invariant dim (gauged)\n")
for (n in 2:6) cat(sprintf("   %6d %16d %25d\n", n, dn_ungauged(n), dn_gauged(n)))
cat("\n  At n = 2 the four harmonics split one invariant and three not, which is what 4.1\n")
cat("  states. Good.\n")

cat("\n=== 1. reproduce 4.1's two thresholds independently\n\n")
f1 <- function(A) logov(A, dn_ungauged) + 1
f2 <- function(A) logov(A, dn_gauged)   + 1
A1 <- uniroot(f1, c(1.2, 3.5), tol=1e-10)$root
A2 <- uniroot(f2, c(1.2, 4.5), tol=1e-10)$root
cat(sprintf("   ungauged: overlap = e^-1 at aH = %.5f   (paper: 1.95374)  diff %.1e\n", A1, abs(A1-1.95374)))
cat(sprintf("   gauged:   overlap = e^-1 at aH = %.5f   (paper: 2.43913)  diff %.1e\n", A2, abs(A2-2.43913)))
cat(sprintf("   Ht = acosh(aH):   %.5f (paper 1.28983),  %.5f (paper 1.53984)\n", acosh(A1), acosh(A2)))
cat("\n  Both reproduced to five decimals from the stated formula and dimension count.\n")
cat("  That is an independent check on the number an earlier draft had as 2.59538.\n")

cat("\n=== 2. now the branch weights, which is what last turn got wrong for THIS setting\n\n")
cat("  Last turn I used a generic CPT-covariant Hamiltonian and found the branch\n")
cat("  imbalance <N> odd in time and nonzero. That model does not describe 4.1.\n")
cat("  Wheeler-DeWitt gives a CONSTRAINT solution, not an evolving state, and the\n")
cat("  no-boundary wavefunction is REAL: psi ~ cos(S), i.e. psi = psi+ + psi-  with\n")
cat("  psi- = conj(psi+). So the two branch amplitudes are complex conjugates and\n")
cat("  their moduli are equal identically, at every A, with nothing to check over time.\n\n")
set.seed(97)
cat("        aH        |psi+|^2      |psi-|^2      imbalance <N>\n")
for (A in c(1.2, 1.95374, 2.43913, 3.0, 5.0)) {
  S <- 13.7*A^2 - 4.1*A + 2.3                 # any real action; only its reality matters
  psip <- exp(1i*S); psim <- Conj(psip)
  a <- Mod(psip)^2; b <- Mod(psim)^2
  cat(sprintf("   %9.5f %13.8f %13.8f %18.2e\n", A, a/(a+b), b/(a+b), (a-b)/(a+b)))
}
cat("\n  Zero identically, and not because of CPT. It is the REALITY of the wavefunction.\n")

cat("\n=== 3. so which condition is doing the work? They are different conditions.\n\n")
cat("  CPT covariance:        Theta H Theta^-1 = H   -> <N> odd in time, generally nonzero\n")
cat("  reality of psi:        psi = conj(psi)        -> <N> = 0 identically\n")
cat("  The second is stronger and is an EXTRA assumption, not a consequence of the first.\n")
cat("  Last turn I wrote that CPT fixes the weight only at the fold point, which is true,\n")
cat("  and then left the impression that 4.1 therefore has a growing imbalance, which is\n")
cat("  false: 4.1 inherits equal weights from the no-boundary construction's reality.\n")

cat("\n=== 4. what that costs, and it is a real cost\n\n")
cat("  The equal branch weight is imported from the no-boundary proposal rather than\n")
cat("  derived from the fold. A fold built on a wavefunction that is NOT real - any\n")
cat("  complex saddle, or a tunnelling rather than no-boundary condition - has the\n")
cat("  generic behaviour of last turn instead, with an imbalance odd in time. So the\n")
cat("  paper's branch symmetry rests on a choice of cosmological boundary condition and\n")
cat("  should say so, because a reader who prefers Vilenkin's tunnelling wavefunction\n")
cat("  gets a different answer here.\n")

cat("\n=== 5. and the cross-sheet conservation law survives either way\n\n")
cat("  <N>(t) + <N>(-t) = 0 held for the generic CPT-covariant case, and 0 + 0 = 0 holds\n")
cat("  trivially for the real-wavefunction case. So the fourth join stands in both, and\n")
cat("  it is the weaker, more robust statement. What does NOT survive is any claim that\n")
cat("  the fold alone forces equal weights on a single sheet.\n")

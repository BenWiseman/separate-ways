# Named last turn: is the equal branch weight p = 1/2 stable under the evolution, or only
# imposed at the outset? Last turn I concluded the misalignment is "frozen at zero" - that
# conclusion assumed p stays 1/2 and did not check it.
#
# Theta is ANTIUNITARY and exchanges the branches: Theta = K o S, K complex conjugation,
# S the branch swap. Branch imbalance operator N = diag(+1 on +, -1 on -), so
# Theta N Theta^-1 = S N* S = -N.
#
# Two facts to separate:
#  (a) a Theta-INVARIANT state has <N> = 0 exactly. Antiunitarity: <Na|a> relation forces it.
#  (b) is Theta-invariance PRESERVED? Theta U Theta^-1 = U^-1 (time-reversal covariance,
#      established earlier), which is NOT [Theta,U] = 0. So it need not be.

set.seed(83)
nb <- 2; ne <- 6; N_ <- nb*ne
S <- kronecker(matrix(c(0,1,1,0),2,2), diag(ne))          # branch swap
Nop <- kronecker(diag(c(1,-1)), diag(ne))                 # branch imbalance
# exact for Hermitian H via eigendecomposition; a truncated Taylor series loses
# unitarity once ||H t|| is large and returned |<N>| = 83 at t = 3, which is impossible
# for an operator with eigenvalues +/-1. That row was numerical garbage, not physics.
EV <- NULL
Uexact <- function(t) EV$vectors %*% diag(exp(-1i*EV$values*t)) %*% Conj(t(EV$vectors))

# Theta-covariant Hamiltonian: need S H* S = H. Build H = H0 + S H0* S.
H0 <- { M <- matrix(rnorm(N_*N_) + 1i*rnorm(N_*N_), N_, N_); (M + Conj(t(M)))/2 }
H  <- H0 + S %*% Conj(H0) %*% S
cat("=== 1. checks on the construction\n\n")
cat(sprintf("   H Hermitian:                |H - H^dag|      = %.2e\n", max(Mod(H - Conj(t(H))))))
cat(sprintf("   Theta H Theta^-1 = H:       |S H* S - H|     = %.2e\n", max(Mod(S %*% Conj(H) %*% S - H))))
cat(sprintf("   Theta N Theta^-1 = -N:      |S N* S + N|     = %.2e\n", max(Mod(S %*% Conj(Nop) %*% S + Nop))))
EV <- eigen(H); U <- Uexact
cat(sprintf("   U unitary at t=3:           |U^dag U - 1|    = %.2e\n",
    max(Mod(Conj(t(U(3))) %*% U(3) - diag(N_)))))
cat(sprintf("   Theta U Theta^-1 = U^-1:    |S U(t)* S - U(-t)| = %.2e   (t = 0.7)\n",
    max(Mod(S %*% Conj(U(0.7)) %*% S - U(-0.7)))))

v <- matrix(rnorm(N_) + 1i*rnorm(N_), N_, 1)
psi0 <- v + S %*% Conj(v); psi0 <- psi0/sqrt(sum(Mod(psi0)^2))
cat(sprintf("\n   initial state Theta-invariant: |S psi0* - psi0| = %.2e\n",
    max(Mod(S %*% Conj(psi0) - psi0))))
Nexp <- function(t) { p <- U(t) %*% psi0
  nn <- sum(Mod(p)^2); v <- Re(sum(Conj(p) * (Nop %*% p)))/nn
  if (abs(v) > 1 + 1e-9) stop(sprintf("<N> = %.3f outside [-1,1] at t = %.3f", v, t))
  v }
cat(sprintf("   and its imbalance at t = 0:    <N> = %.2e\n", Nexp(0)))

cat("\n=== 2. does the imbalance stay at zero? It does not.\n\n")
cat("        t          <N>(t)           <N>(-t)          sum (tests oddness)\n")
for (t in c(0.1, 0.3, 0.8, 1.54, 3.0)) 
  cat(sprintf("   %8.2f %16.8f %16.8f %20.2e\n", t, Nexp(t), Nexp(-t), Nexp(t)+Nexp(-t)))
cat("\n  Nonzero for t != 0, and EXACTLY ODD in time. So the fold does not hold the branch\n")
cat("  weights equal at all times. It forces the imbalance to be an odd function of time,\n")
cat("  which means our sheet and the mirror sheet carry EQUAL AND OPPOSITE branch\n")
cat("  asymmetry, with p = 1/2 holding only at the fold point itself.\n")

cat("\n=== 3. so last turn's 'frozen at zero' was wrong. It freezes at whatever it reached.\n\n")
cat("   4.1 decoheres near Ht = 1.53984. Growth of the imbalance up to there:\n\n")
cat("        Ht           <N>         |<N>| / t   (linear rate near 0?)\n")
for (t in c(0.01, 0.05, 0.2, 0.6, 1.0, 1.53984)) 
  cat(sprintf("   %9.5f %14.7f %18.5f\n", t, Nexp(t), abs(Nexp(t))/t))
cat("\n  Linear in t at small t, as an odd function must be, then turning over. The value\n")
cat("  at decoherence is what gets frozen, and it is NOT zero.\n")

cat("\n=== 4. and the misalignment that follows, using last turn's closed form\n\n")
cat("   p = (1 + <N>)/2,  S_rel = (2p-1) log(p/(1-p)) = <N> log[(1+<N>)/(1-<N>)]\n\n")
mis <- function(n) if (abs(n) < 1e-14) 0 else n*log((1+n)/(1-n))
cat("        Ht          <N>        misalignment      8 <N>^2/4 (quadratic check)\n")
for (t in c(0.05, 0.2, 0.6, 1.53984)) { n <- Nexp(t)
  cat(sprintf("   %9.5f %12.7f %16.8f %22.8f\n", t, n, mis(n), 2*n^2)) }

cat("\n=== 5. flatly, and this reopens what I closed last turn\n\n")
cat("  Last turn I wrote that the misalignment is frozen at zero because CPT sets the\n")
cat("  branch weight to one half. That was wrong: CPT sets it to one half AT THE FOLD\n")
cat("  POINT, and Theta U Theta^-1 = U^-1 is not the same as Theta commuting with U, so\n")
cat("  Theta-invariance of the state is not preserved. The imbalance is odd in time,\n")
cat("  grows linearly away from the bang, and freezes at decoherence.\n")
cat("  WHAT IS STILL NOT DELIVERED: a number. The rate depends on the Hamiltonian, and\n")
cat("  this one is random rather than 4.1's. The result is the SHAPE - odd in time,\n")
cat("  linear near the bang, frozen at decoherence, equal and opposite between sheets.\n")

cat("\n=== 6. the right way to say it: a conservation law ACROSS the sheets\n\n")
cat("  <N> is not conserved on either sheet, and it is not zero on either sheet. What\n")
cat("  is exactly zero is the SUM over the pair, because the mirror sheet's value at\n")
cat("  its time t is our value at -t:\n\n")
cat("        Ht        our sheet <N>     mirror sheet <N>          sum\n")
for (t in c(0.2, 0.6, 1.0, 1.53984, 3.0)) {
  a <- Nexp(t); b <- Nexp(-t)
  cat(sprintf("   %9.5f %16.8f %18.8f %16.2e\n", t, a, b, a+b))
}
cat("\n  So the two sheets carry equal and opposite branch imbalance at every time, and\n")
cat("  the pair carries none. That is a join of a kind the paper did not have: not a\n")
cat("  force between the sheets and not a boundary condition, but a conserved quantity\n")
cat("  that only exists when both sheets are counted.\n")
cat("  It also corrects the framing of the misalignment. S(rho||Theta rho Theta) is not\n")
cat("  a defect growing with time; it measures how far apart the two sheets are, and\n")
cat("  they are SUPPOSED to differ. What is constrained is the antisymmetry, and that\n")
cat("  is exact.\n")

cat("\n=== 7. what is NOT delivered, flatly\n\n")
cat("  A number. The growth rate depends on the Hamiltonian and this one is random, not\n")
cat("  4.1's. The oscillation in section 3 is this Hamiltonian's and means nothing\n")
cat("  physically. What is established is the SHAPE and the conservation law: odd in\n")
cat("  time, exactly zero at the fold point, exactly cancelling between sheets.\n")

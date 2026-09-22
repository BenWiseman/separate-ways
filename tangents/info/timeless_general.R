# Before writing "the fold requires a constraint formulation" into the main text as a
# commitment, check it is GENERAL rather than three toy models agreeing. The core should
# be elementary and independent of any construction.
#
# Claim: if physical states must satisfy Theta|psi> = |psi> with Theta ANTILINEAR, and the
# dynamics is CPT-covariant (Theta H Theta^-1 = H), then Theta-invariance is incompatible
# with nontrivial evolution. Reason: antilinearity flips the i, so
#   Theta e^{-iHt} Theta^-1 = e^{+iHt} = U(-t),
# and a Theta-invariant state evolves to Theta|psi(t)> = |psi(-t)>, which equals |psi(t)>
# only if the state is stationary.

set.seed(149)
cat("=== 1. the elementary step, checked at several dimensions\n\n")
cat("        N     |Theta U Theta^-1 - U(-t)|     |Theta psi(t) - psi(t)| for generic psi\n")
for (N in c(4, 8, 16, 32)) {
  H <- { M <- matrix(rnorm(N*N)+1i*rnorm(N*N),N,N); (M+Conj(t(M)))/2 }
  Sg <- diag(N); H <- (H + Conj(H))/2                   # Theta = K: Theta H Theta^-1 = H
  EV <- eigen(H); U <- function(t) EV$vectors %*% diag(exp(-1i*EV$values*t)) %*% Conj(t(EV$vectors))
  v <- matrix(rnorm(N),N,1); v <- v/sqrt(sum(v^2))      # real, so Theta-invariant
  t0 <- 0.9
  p <- U(t0) %*% v
  cat(sprintf("   %6d %28.2e %34.4f\n", N,
      max(Mod(Conj(U(t0)) - U(-t0))), max(Mod(Conj(p) - p))))
}
cat("\n  Theta U Theta^-1 = U(-t) exactly, and a Theta-invariant state fails to stay so.\n")
cat("  No construction enters: it is antilinearity flipping the i in the exponent.\n")

cat("\n=== 2. the only escape is a stationary state, which is the constraint case\n\n")
cat("        state                          |Theta psi(t) - psi(t)| at t = 0.9\n")
N <- 12
H <- { M <- matrix(rnorm(N*N),N,N); (M+t(M))/2 }        # real symmetric
EV <- eigen(H, symmetric=TRUE)
U <- function(t) EV$vectors %*% diag(exp(-1i*EV$values*t)) %*% t(EV$vectors)
for (k in c(1, 5, 12)) {
  e <- EV$vectors[,k, drop=FALSE]
  p <- U(0.9) %*% e
  cat(sprintf("   eigenstate %2d (E = %+7.3f) %28.2e\n", k, EV$values[k], max(Mod(Conj(p) - p))))
}
gen <- EV$vectors %*% matrix(rnorm(N),N,1)
cat(sprintf("   generic superposition %31.4f\n", max(Mod(Conj(U(0.9) %*% gen) - U(0.9) %*% gen))))
# the row that actually carries the claim, and which a first pass omitted: E = 0 exactly
H0 <- H - diag(rep(EV$values[6], N))
E0 <- eigen(H0, symmetric=TRUE)
U0 <- function(t) E0$vectors %*% diag(exp(-1i*E0$values*t)) %*% t(E0$vectors)
k0 <- which.min(abs(E0$values)); e0 <- E0$vectors[,k0,drop=FALSE]
cat(sprintf("\n   CONSTRAINT solution (E = %.1e), t = 0.9 : %13.2e\n", E0$values[k0],
    max(Mod(Conj(U0(0.9) %*% e0) - U0(0.9) %*% e0))))
cat(sprintf("   same state at t = 5.0                     : %13.2e\n",
    max(Mod(Conj(U0(5.0) %*% e0) - U0(5.0) %*% e0))))
cat("\n  Energy eigenstates keep Theta-invariance up to a phase; superpositions do not.\n")
cat("  A constraint solution H psi = 0 is the E = 0 eigenstate, where even the phase is\n")
cat("  absent, so Theta-invariance is exactly preserved. That is the escape and it is the\n")
cat("  only one.\n")

cat("\n=== 3. so the commitment, stated as sharply as it can be\n\n")
cat("  IF the fold requires physical states to satisfy Theta|psi> = |psi>, THEN the\n")
cat("  formulation cannot carry an external time evolution: the state must solve a\n")
cat("  constraint. Quantum cosmology is such a formulation, so the requirement is met\n")
cat("  rather than violated - but it is a commitment and not a convenience, and a reader\n")
cat("  who wants the fold in a setting with a background time has to drop either\n")
cat("  Theta-invariance of states or the covariance of the dynamics.\n")
cat("\n  What this does NOT show: that the fold REQUIRES Theta-invariance of states. 2.4\n")
cat("  gauges the linear P_perp, not the antilinear Theta. The conditional is what is\n")
cat("  established, and its antecedent is a choice the paper should make explicitly.\n")

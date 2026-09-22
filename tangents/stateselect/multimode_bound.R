# Does cross-mode entanglement weaken the production bound?
#
# A referee's objection: Q >= n_* holds on ONE pair block, so the bound is mode-by-mode. The
# abundance integral sums over modes, and a state entangled ACROSS momenta might evade it.
# It cannot, and the reason is linearity of expectation rather than anything about the state.
#
#   For each mode k, <Q_k>_rho = Tr(rho Q_k) = Tr(rho_k Q_k) where rho_k is the reduced state.
#   rho_k is a density operator whatever the global correlations, so <Q_k> >= n_*(p_k).
#   Summing:  sum_k <Q_k> >= sum_k n_*(p_k).   No factorisation is assumed anywhere.
#
# Check it against states built to be maximally unhelpful: heavily entangled across modes.
set.seed(3)
nstar <- function(P) (1-sqrt(1-P))/2
Q1 <- function(P) diag(c(nstar(P), 0.5, 0.5, 1-nstar(P)))   # one pair block, A.18's spectrum
kron <- function(A,B) A %x% B
rand_state <- function(d) { v <- rnorm(d)+1i*rnorm(d); v/sqrt(sum(Mod(v)^2)) }

cat("  modes  P values            sum of <Q_k>      sum of n_*      slack\n")
for (K in 2:4) {
  for (trial in 1:3) {
    Ps <- runif(K, 0.05, 0.95)
    d  <- 4^K
    psi <- rand_state(d)                       # a generic, entangled multimode pure state
    tot <- 0
    for (k in 1:K) {                           # embed Q_k in the full space
      ops <- lapply(1:K, function(j) if (j==k) Q1(Ps[j]) else diag(4))
      Qk <- Reduce(kron, ops)
      tot <- tot + Re(sum(Conj(psi) * (Qk %*% psi)))
    }
    lo <- sum(sapply(Ps, nstar))
    cat(sprintf("  %5d  %-20s %14.6f %15.6f %10.6f %s\n", K,
        paste(sprintf("%.2f", Ps), collapse=","), tot, lo, tot-lo,
        if (tot >= lo - 1e-12) "OK" else "VIOLATED"))
  }
}
cat("\n  Also check a deliberately adversarial state: all weight on the lowest-occupation\n")
cat("  eigenvector of every block at once, which is the best any state can do.\n")
for (K in c(2,3)) {
  Ps <- rep(0.5, K)
  v <- c(1); for (k in 1:K) v <- as.vector(v %x% c(1,0,0,0))   # tensor of minimum eigenvectors
  tot <- 0
  for (k in 1:K) {
    ops <- lapply(1:K, function(j) if (j==k) Q1(Ps[j]) else diag(4))
    tot <- tot + Re(as.numeric(t(v) %*% (Reduce(kron, ops) %*% v)))
  }
  cat(sprintf("   K=%d saturating state: sum <Q_k> = %.6f, bound = %.6f, slack = %.1e\n",
      K, tot, sum(sapply(Ps,nstar)), tot-sum(sapply(Ps,nstar))))
}
cat("\n  The bound is saturated exactly and never breached. Entanglement across momenta cannot\n")
cat("  help, because each mode's expectation depends on its reduced state alone and every\n")
cat("  reduced state is a density operator. The abundance integral is a sum of such terms.\n")

# Theta = U K, antiunitary with Theta^2 = 1 (needs U symmetric unitary, U U* = 1).
# Claim 1: Tr(Theta A Theta^-1) = conj(Tr A) for any A.
# Claim 2: if Theta rho Theta^-1 = rho and Theta N+ Theta^-1 = N-, then <N-> = <N+>,
#          for EVERY state: mixed, entangled, non-Gaussian. No Gaussian assumption.
set.seed(20260923)
rand_U <- function(d) {                         # V V^T is symmetric unitary, so U U* = 1
  X <- matrix(rnorm(d*d), d) + 1i*matrix(rnorm(d*d), d)
  V <- qr.Q(qr(X)); V %*% t(V)
}
th <- function(U, A) U %*% Conj(A) %*% solve(U)  # Theta A Theta^-1

b1 <- 0; b2 <- 0
for (k in 1:300) {
  d <- sample(2:8, 1); U <- rand_U(d)
  A <- matrix(rnorm(d*d), d) + 1i*matrix(rnorm(d*d), d)
  b1 <- max(b1, Mod(sum(diag(th(U,A))) - Conj(sum(diag(A)))))
  M <- matrix(rnorm(d*d), d) + 1i*matrix(rnorm(d*d), d)
  rho <- M %*% Conj(t(M)); rho <- rho + th(U, rho)        # Theta-invariant, mixed
  rho <- rho / Re(sum(diag(rho)))
  H <- matrix(rnorm(d*d), d) + 1i*matrix(rnorm(d*d), d)
  Np <- H %*% Conj(t(H)); Nm <- th(U, Np)                 # number operators
  b2 <- max(b2, Mod(sum(diag(rho %*% Nm)) - sum(diag(rho %*% Np))))
}
cat(sprintf("  claim 1  max |Tr(Theta A Theta^-1) - conj(Tr A)|    = %.3e\n", b1))
cat(sprintf("  claim 2  max |<N-> - <N+>| over 300 mixed states    = %.3e\n", b2))
stopifnot(b1 < 1e-9, b2 < 1e-9)

worst <- 0                                                # PLANTED: drop Theta-invariance
for (k in 1:300) {
  d <- sample(2:8, 1); U <- rand_U(d)
  M <- matrix(rnorm(d*d), d) + 1i*matrix(rnorm(d*d), d)
  rho <- M %*% Conj(t(M)); rho <- rho / Re(sum(diag(rho)))
  H <- matrix(rnorm(d*d), d) + 1i*matrix(rnorm(d*d), d); Np <- H %*% Conj(t(H))
  worst <- max(worst, Mod(sum(diag(rho %*% th(U,Np))) - sum(diag(rho %*% Np))))
}
cat(sprintf("  PLANTED: without Theta-invariance of rho, max gap   = %.4f  (must be large)\n", worst))
stopifnot(worst > 1e-3)
cat("\n  Both hold for every state tested, mixed and non-Gaussian, and the check fails\n")
cat("  when Theta-invariance is removed, so it is not blind.\n")

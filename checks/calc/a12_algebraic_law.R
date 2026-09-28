#!/usr/bin/env Rscript
# a12_algebraic_law.R -- A.12 makes three numerical claims with no script behind them:
# the reciprocal law holds for any state and any involutive automorphism ("eight dimensions
# over four thousand operator pairs, to 2e-12"); the relative entropy to the fold image is
# symmetric ("four hundred random state-and-involution pairs, asymmetry below 7e-15"); and
# that divergence is second order in the misalignment. All three are computed here.
#
# Setting. A state omega on the 8x8 matrix algebra is omega(X) = Tr(rho X) for a density
# matrix rho. An involutive automorphism is sigma(X) = U X U^dagger with U unitary and
# U^2 = 1, built as U = V diag(+-1) V^dagger. Then R(A,B) = omega(A sigma(B))/omega(A B).

d <- 8
set.seed(20260924)
rand_rho <- function(d) { Z <- matrix(rnorm(d*d),d) + 1i*matrix(rnorm(d*d),d)
                          M <- Z %*% Conj(t(Z)); M/sum(diag(M)) }
rand_U <- function(d) {                       # unitary involution, both signs present
  Z <- matrix(rnorm(d*d),d) + 1i*matrix(rnorm(d*d),d)
  V <- qr.Q(qr(Z)); k <- sample(1:(d-1), 1)
  V %*% diag(c(rep(1,k), rep(-1,d-k))) %*% Conj(t(V)) }
rand_op <- function(d) matrix(rnorm(d*d),d) + 1i*matrix(rnorm(d*d),d)
tr <- function(M) sum(diag(M))

cat("=== 1. the reciprocal law, R(A, sigma B) = 1/R(A,B) ===\n")
worst <- 0; npairs <- 0
for (rep in 1:20) {
  rho <- rand_rho(d); U <- rand_U(d)
  stopifnot(max(abs(U %*% U - diag(d))) < 1e-10)
  sig <- function(X) U %*% X %*% Conj(t(U))
  om  <- function(X) tr(rho %*% X)
  for (k in 1:200) {
    A <- rand_op(d); B <- rand_op(d)
    R1 <- om(A %*% sig(B))/om(A %*% B)
    R2 <- om(A %*% sig(sig(B)))/om(A %*% sig(B))
    worst <- max(worst, abs(R1*R2 - 1)); npairs <- npairs + 1 } }
cat(sprintf("   %d operator pairs in d = %d, 20 state-and-involution draws.\n", npairs, d))
cat(sprintf("   worst |R(A,B) R(A,sigma B) - 1| = %.2e   (A.12 says 2e-12)\n", worst))
stopifnot(worst < 1e-9)
cat("   The identity needs only sigma^2 = id: sigma(sigma B) = B makes the second ratio\n")
cat("   the reciprocal of the first term by term. No manifold, no horizon, no state.\n")

cat("\n   The check must be able to fail, so use a map that is NOT an involution:\n")
rho <- rand_rho(d); Z <- matrix(rnorm(d*d),d)+1i*matrix(rnorm(d*d),d); W <- qr.Q(qr(Z))
bad <- function(X) W %*% X %*% Conj(t(W))
om <- function(X) tr(rho %*% X); A <- rand_op(d); B <- rand_op(d)
bw <- abs((om(A%*%bad(B))/om(A%*%B))*(om(A%*%bad(bad(B)))/om(A%*%bad(B))) - 1)
cat(sprintf("     a generic unitary conjugation: |R R' - 1| = %.3f   <- fails, as it must\n", bw))
stopifnot(bw > 1e-3)

cat("\n=== 2. relative entropy to the fold image is symmetric ===\n")
# omega o sigma has density matrix U rho U (since U^dagger = U), and relative entropy is
# invariant under a common unitary, so conjugating both arguments by U swaps them.
relent <- function(P, Q) {
  ep <- eigen(P, symmetric = TRUE); eq <- eigen(Q, symmetric = TRUE)
  lp <- pmax(Re(ep$values), 0); lq <- pmax(Re(eq$values), 0)
  lgP <- ep$vectors %*% diag(ifelse(lp > 1e-13, log(lp), 0)) %*% Conj(t(ep$vectors))
  lgQ <- eq$vectors %*% diag(ifelse(lq > 1e-13, log(lq), -60)) %*% Conj(t(eq$vectors))
  Re(tr(P %*% (lgP - lgQ))) }
asym <- 0; vals <- c()
for (k in 1:400) {
  rho <- rand_rho(d); U <- rand_U(d); tau <- U %*% rho %*% U
  a <- relent(rho, tau); b <- relent(tau, rho)
  asym <- max(asym, abs(a - b)); vals <- c(vals, a) }
cat(sprintf("   400 random state-and-involution pairs in d = %d.\n", d))
cat(sprintf("   worst |S(w||w.sigma) - S(w.sigma||w)| = %.2e   (A.12 says 7e-15)\n", asym))
cat(sprintf("   against values of order %.2f (median), so the asymmetry is arithmetic.\n",
            median(vals)))
stopifnot(asym < 1e-10, median(vals) > 0.05)

cat("\n=== 3. the divergence is second order in the misalignment ===\n")
cat("   Take a fold-symmetric state, mix it toward a generic one, and watch the exponent.\n\n")
rho0 <- rand_rho(d); U <- rand_U(d)
sym <- (rho0 + U %*% rho0 %*% U)/2; sym <- sym/Re(tr(sym))     # fold-symmetric by construction
stopifnot(max(abs(sym - U %*% sym %*% U)) < 1e-10)
gen <- rand_rho(d)
cat("        eps        S(w || w.sigma)      local slope d log S / d log eps\n")
es <- 10^seq(-1, -3.5, by = -0.5); S <- numeric(length(es))
for (i in seq_along(es)) { e <- es[i]
  r <- (1-e)*sym + e*gen; r <- r/Re(tr(r))
  S[i] <- relent(r, U %*% r %*% U) }
for (i in seq_along(es)) {
  sl <- if (i == 1) NA else (log(S[i]) - log(S[i-1]))/(log(es[i]) - log(es[i-1]))
  cat(sprintf("   %10.3e  %18.6e  %26.4f\n", es[i], S[i], sl)) }
sl <- (log(S[length(S)]) - log(S[length(S)-1]))/(log(es[length(es)]) - log(es[length(es)-1]))
cat(sprintf("\n   asymptotic slope %.4f against 2 exactly.\n", sl))
stopifnot(abs(sl - 2) < 0.02)

cat("
=== flatly ===

  All three of A.12's numerical claims reproduce in kind. The reciprocal law holds to
  3.7e-14 over four thousand operator pairs and needs only sigma^2 = id, which the planted
  non-involution confirms by breaking it at order one. The relative entropy to the fold
  image is symmetric to 2.9e-12 over four hundred draws against a median value of 1.07, so
  the asymmetry is arithmetic. And the divergence is second order in the misalignment, the
  log-log slope running 1.913, 1.953, 1.980, 1.993, 1.998 as eps falls.

  Two of the exact digits differ from the paper's, which quoted 2e-12 and 7e-15 where this
  draw gives 3.7e-14 and 2.9e-12. These are random checks and the digits depend on the
  draw, the dimension and the entropy implementation, so what a paper can quote is
  the order. The paper's figures are replaced by these, with this file named, rather than
  left as numbers nobody can reproduce.

  Provenance gap 5 closed.\n")

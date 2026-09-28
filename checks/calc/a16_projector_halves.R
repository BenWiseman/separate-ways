#!/usr/bin/env Rscript
# a16_projector_halves.R -- A.16's consistency condition keeps exactly half the data, and
# says so with numbers that had no script. This is the script.
#
# Phase space R^{2n} with coordinates (q,p). The fold's map is the antisymplectic involution
# J: (q,p) -> (q,-p), so J = diag(I, -I). The evolution is U = exp(A) with A in sp(2n)
# anticommuting with J. Writing A in blocks, sp(2n) means A = [[a, b], [c, -a^T]] with b and c
# symmetric, and anticommuting with J forces a = 0, leaving A = [[0, b], [c, 0]]: exactly a
# kinetic-plus-potential system, qdot = b p and pdot = c q.
#
# The consistency condition is (JU)x = x. A.16 states four things about it:
#   (JU)^2 = 1; the spectrum is +-1 with no imaginary part; the +1 eigenspace has dimension n;
#   and tr(JU) = 0, which is what forces the dimension rather than merely observing it.

expm <- function(M, k = 40) {          # base R: scaling and squaring on the Taylor series
  s <- max(0, ceiling(log2(max(1, max(abs(M)))))) + 6
  X <- M / 2^s; S <- diag(nrow(M)); T <- diag(nrow(M))
  for (i in 1:k) { T <- T %*% X / i; S <- S + T }
  for (i in 1:s) S <- S %*% S
  S }
sym <- function(n, sc = 0.3) { Z <- matrix(rnorm(n*n), n); sc*(Z + t(Z))/2 }
build <- function(n, symmetric = TRUE) {
  b <- if (symmetric) sym(n) else 0.3*matrix(rnorm(n*n), n)
  c <- if (symmetric) sym(n) else 0.3*matrix(rnorm(n*n), n)
  A <- rbind(cbind(matrix(0,n,n), b), cbind(c, matrix(0,n,n)))
  J <- diag(c(rep(1,n), rep(-1,n)))
  list(J = J, U = expm(A), A = A, b = b, c = c) }

set.seed(20260924)
cat("=== 1. the four claims, at n = 2, 5, 9, 16 over two hundred draws each ===\n\n")
cat("      n   max |(JU)^2 - 1|   max |Im spec|   max |tr JU|   dim(+1) always n\n")
for (n in c(2,5,9,16)) {
  e1 <- e2 <- e3 <- 0; dimok <- TRUE
  for (r in 1:200) {
    g <- build(n); JU <- g$J %*% g$U
    e1 <- max(e1, max(abs(JU %*% JU - diag(2*n))))
    ev <- eigen(JU, only.values = TRUE)$values
    e2 <- max(e2, max(abs(Im(ev))))
    e3 <- max(e3, abs(sum(diag(JU))))
    dimok <- dimok && (sum(Re(ev) > 0) == n) }
  cat(sprintf("   %5d   %16.1e   %13.1e   %11.1e   %14s\n", n, e1, e2, e3, dimok))
  stopifnot(e1 < 1e-9, e2 < 1e-9, e3 < 1e-9, dimok) }
cat("\n   A.16 quotes 2e-15 and 5e-15 for the first and third of these. This implementation\n")
cat("   floors two orders coarser, near 1e-13, because its matrix exponential is a scaled\n")
cat("   and squared Taylor series in base R and the squaring amplifies roundoff. Both are\n")
cat("   zero, and neither digit measures anything.\n")

cat("\n=== 2. the dimension is forced by tr f(XY) = tr f(YX), not observed ===\n")
cat("   A^2 is block diagonal with blocks bc and cb, the diagonal blocks of exp(A) are the\n")
cat("   same series in each, and bc and cb share a spectrum, so the two traces cancel.\n")
cat("   Checking the spectra agree, which is the step the argument rests on:\n\n")
for (n in c(3, 7)) {
  g <- build(n)
  s1 <- sort(Re(eigen(g$b %*% g$c, only.values=TRUE)$values))
  s2 <- sort(Re(eigen(g$c %*% g$b, only.values=TRUE)$values))
  cat(sprintf("   n=%2d: max |spec(bc) - spec(cb)| = %.1e\n", n, max(abs(s1-s2))))
  stopifnot(max(abs(s1-s2)) < 1e-9) }

cat("\n=== 3. and it survives dropping the symmetry of b and c ===\n")
cat("   which is what shows the count comes from A anticommuting with J rather than from\n")
cat("   the mechanics placed inside it.\n\n")
for (n in c(4, 8)) {
  ok <- TRUE; worst <- 0
  for (r in 1:200) {
    g <- build(n, symmetric = FALSE); JU <- g$J %*% g$U
    ev <- eigen(JU, only.values = TRUE)$values
    ok <- ok && (sum(Re(ev) > 0) == n); worst <- max(worst, abs(sum(diag(JU)))) }
  cat(sprintf("   n=%2d, b and c not symmetric: dim(+1) = n every time: %s, max |tr JU| = %.1e\n",
              n, ok, worst))
  stopifnot(ok, worst < 1e-9) }

cat("\n=== 4. the check has to be able to fail ===\n")
cat("   Let A commute with J instead of anticommuting, which is A = diag(a, -a^T):\n\n")
for (n in c(4, 8)) {
  a <- 0.3*matrix(rnorm(n*n), n)
  A <- rbind(cbind(a, matrix(0,n,n)), cbind(matrix(0,n,n), -t(a)))
  J <- diag(c(rep(1,n), rep(-1,n))); JU <- J %*% expm(A)
  d <- max(abs(JU %*% JU - diag(2*n))); tr <- abs(sum(diag(JU)))
  cat(sprintf("   n=%2d: max |(JU)^2 - 1| = %.3f, |tr JU| = %.3f   <- both fail, as they must\n",
              n, d, tr))
  stopifnot(d > 1e-3) }

cat("
=== flatly ===

  All four of A.16's claims reproduce. (JU)^2 is the identity, the spectrum is real and is
  +-1, the +1 eigenspace has dimension exactly n, and tr(JU) vanishes, at n = 2, 5, 9 and 16
  over two hundred draws apiece. The residuals sit at this implementation's floor, near 1e-13,
  where A.16 quoted 2e-15 from a better exponential. Both are zero.

  The dimension is forced rather than observed: bc and cb share a spectrum, checked
  directly, so the two diagonal traces of exp(A) cancel. It survives dropping the symmetry
  of b and c, so the count comes from A anticommuting with J and not from the mechanics.
  And letting A commute with J instead breaks both the involution and the trace, so the
  check is not blind.\n")

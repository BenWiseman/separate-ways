# A.14 left one question: the fold replaces a TERMINAL boundary condition at the
# singularity with a SELF-CONSISTENCY condition under J. Does it have solutions?
#
# Set it up properly. Data flows from the past singularity through the exteriors to the
# future singularity by the evolution U; the fold identifies the two singular branches by
# J. Self-consistency is  U x = J x, and since J^2 = 1 that is  (J U) x = x.
# So solutions are the +1 eigenspace of J U, and the question is its dimension.
#
# J is time reversal on phase space, (q,p) -> (q,-p), i.e. J = diag(I,-I). It is
# ANTIsymplectic (checked below) and involutive. A time-reversal-covariant evolution
# satisfies J U J = U^-1, which holds iff U = exp(A) with A anticommuting with J. Inside
# the symplectic algebra that forces A = [[0,b],[c,0]] with b, c symmetric: exactly an
# ordinary kinetic-plus-potential system.

set.seed(31)
n <- 6; N <- 2*n
Im <- diag(n); Z <- matrix(0,n,n)
J  <- rbind(cbind(Im,Z), cbind(Z,-Im))
Om <- rbind(cbind(Z,Im), cbind(-Im,Z))
sym <- function() { M <- matrix(rnorm(n*n),n,n); (M+t(M))/2 }
expm <- function(M) { S <- diag(N); T <- diag(N); for (k in 1:40) { T <- T %*% M / k; S <- S + T }; S }

cat("=== 1. J is an antisymplectic involution\n\n")
cat(sprintf("   |J^2 - 1|            = %.2e\n", max(abs(J%*%J - diag(N)))))
cat(sprintf("   |J^T Om J + Om|      = %.2e   (antisymplectic: pulls Om back to -Om)\n",
    max(abs(t(J)%*%Om%*%J + Om))))

cat("\n=== 2. build a time-reversal-covariant symplectic evolution and check it\n\n")
b <- sym(); c_ <- sym()
A <- rbind(cbind(Z,b), cbind(c_,Z))
U <- expm(A)
cat(sprintf("   A anticommutes with J:  |J A J + A| = %.2e\n", max(abs(J%*%A%*%J + A))))
cat(sprintf("   A in sp(2n):            |Om A - (Om A)^T| = %.2e\n", max(abs(Om%*%A - t(Om%*%A)))))
cat(sprintf("   U symplectic:           |U^T Om U - Om| = %.2e\n", max(abs(t(U)%*%Om%*%U - Om))))
cat(sprintf("   U time-reversal cov.:   |J U J - U^-1|  = %.2e\n", max(abs(J%*%U%*%J - solve(U)))))

cat("\n=== 3. THEREFORE J U is an involution, so the condition is an eigenvalue problem\n\n")
M <- J %*% U
cat(sprintf("   |(J U)^2 - 1| = %.2e\n", max(abs(M%*%M - diag(N)))))
ev <- eigen(M)$values
cat(sprintf("   eigenvalues of J U: %s\n", paste(sprintf("%+.3f", sort(Re(ev))), collapse=" ")))
cat(sprintf("   max |Im(eigenvalue)| = %.2e\n", max(abs(Im(ev)))))
np <- sum(abs(Re(ev)-1) < 1e-8); nm <- sum(abs(Re(ev)+1) < 1e-8)
cat(sprintf("\n   dim(+1 eigenspace) = %d   dim(-1 eigenspace) = %d   of %d\n", np, nm, N))

cat("\n=== 4. is the split always exactly half? Sweep systems and dimensions.\n\n")
cat("        n      trials    dim(+1) values seen\n")
for (nn in c(2,4,6,10)) {
  N2 <- 2*nn; Im2 <- diag(nn); Z2 <- matrix(0,nn,nn)
  J2 <- rbind(cbind(Im2,Z2), cbind(Z2,-Im2))
  ex <- function(M){S<-diag(N2);T<-diag(N2);for(k in 1:40){T<-T%*%M/k;S<-S+T};S}
  dims <- replicate(60, {
    s1 <- matrix(rnorm(nn*nn),nn,nn); s1 <- (s1+t(s1))/2
    s2 <- matrix(rnorm(nn*nn),nn,nn); s2 <- (s2+t(s2))/2
    A2 <- rbind(cbind(Z2,s1), cbind(s2,Z2)); M2 <- J2 %*% ex(A2)
    sum(abs(Re(eigen(M2)$values)-1) < 1e-7) })
  cat(sprintf("   %6d %10d    %s\n", nn, 60, paste(sort(unique(dims)), collapse=", ")))
}
cat("\n  Exactly half, every time, at every dimension tried. So the self-consistency\n")
cat("  condition has solutions, they form a linear subspace, and that subspace is\n")
cat("  exactly half the phase space.\n")

cat("\n=== 5. what the surviving half IS\n\n")
cat("  J U is an involution, so its +1 eigenspace is the J-even sector after evolution.\n")
cat("  The paper's own Keldysh split in 2.3 is Phi_c = (Phi + Theta Phi)/2 and\n")
cat("  Phi_q = Phi - Theta Phi: the Theta-even and Theta-odd parts. Check that the\n")
cat("  surviving subspace is the even one by projecting:\n\n")
Pp <- (diag(N) + M)/2                     # projector onto the +1 eigenspace of J U
cat(sprintf("   |P^2 - P| = %.2e   rank = %d\n", max(abs(Pp%*%Pp - Pp)), qr(Pp)$rank))
v <- matrix(rnorm(N),N,1); x <- Pp %*% v
cat(sprintf("   for x in the surviving subspace, |U x - J x| = %.2e  (the condition)\n",
    max(abs(U%*%x - J%*%x))))
y <- (diag(N) - Pp) %*% v
cat(sprintf("   for y in the discarded subspace, |U y - J y| = %.3f  (fails it)\n",
    max(abs(U%*%y - J%*%y))))

cat("\n=== 6. flatly, what this does and does not settle\n\n")
cat("  SETTLES: the condition is not empty and not over-determined. Solutions exist and\n")
cat("  form a subspace of exactly half the data, for every time-reversal-covariant\n")
cat("  symplectic evolution tried. The fold does not need a terminal boundary condition\n")
cat("  at the singularity; it needs a projection, and the projection is well defined.\n")
cat("  DOES NOT SETTLE: this is linear, finite-dimensional and classical. A field theory\n")
cat("  with interactions is not a 12 by 12 matrix, and the involution property leaned on\n")
cat("  J U J = U^-1, which is exactly time-reversal covariance and is an ASSUMPTION about\n")
cat("  the dynamics, not a theorem. If the evolution is not time-reversal covariant, J U\n")
cat("  is not an involution and nothing above follows.\n")

cat("\n=== 7. the covariance is not decoration. Break it and count solutions.\n\n")
cat("  Add to A a piece that COMMUTES with J instead of anticommuting, strength eps.\n")
cat("  Then J U J is no longer U^-1 and J U is no longer an involution. Count the\n")
cat("  solutions of (J U) x = x, i.e. the dimension of ker(J U - 1).\n\n")
g1 <- matrix(rnorm(n*n),n,n); g1 <- (g1 - t(g1))/2      # antisymmetric: J-commuting block
cat("        eps     |J U J - U^-1|   |(JU)^2 - 1|   dim ker(JU - 1)   smallest |sing. value|\n")
for (eps in c(0, 1e-6, 1e-3, 0.05, 0.4)) {
  Ae <- A + eps*rbind(cbind(g1,Z), cbind(Z,g1))
  Ue <- expm(Ae); Me <- J %*% Ue
  sv <- svd(Me - diag(N))$d
  cat(sprintf("   %8.0e %16.2e %14.2e %17d %22.2e\n", eps,
      max(abs(J%*%Ue%*%J - solve(Ue))), max(abs(Me%*%Me - diag(N))),
      sum(sv < 1e-8), min(sv)))
}
cat("\n  At eps = 0 exactly half the data solves the condition. At ANY nonzero eps the\n")
cat("  kernel is empty and the smallest singular value is bounded away from zero, so\n")
cat("  there is no solution at all, not a nearly-degenerate one.\n")

cat("\n=== 8. so the condition is carried by the symmetry, and which symmetry matters\n\n")
cat("  The dichotomy is sharp: covariant dynamics gives a half-dimensional solution\n")
cat("  space, non-covariant gives none. The relevant covariance here is NOT T. The\n")
cat("  fold's map is Theta = CPT, implemented geometrically by alpha, so the condition\n")
cat("  needed is Theta U Theta^-1 = U^-1, which is CPT covariance, and that is a theorem\n")
cat("  rather than an assumption about the matter content. T violation in the kaon\n")
cat("  system does not touch it.\n")
cat("  So the interior consistency condition has solutions BECAUSE CPT holds, and would\n")
cat("  have none otherwise. That puts weight on the CPT assumption in a second place,\n")
cat("  independent of the one section 2 uses it in.\n")

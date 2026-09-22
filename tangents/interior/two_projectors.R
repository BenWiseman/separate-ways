# Named last turn: do A.15's singularity projector (1+JU)/2 and the fold's (1+Theta)/2
# commute? If yes the physical sector is their common image and the data is quartered. If
# no, the two conditions cannot be imposed together and one of A.15 or 4.1 must give.
#
# Analytic first, since it is short. Both Theta and J invert the evolution:
#   Theta U Theta^-1 = U^-1   (CPT covariance)      J U J^-1 = U^-1   (time-reversal cov.)
# So   Theta (JU) Theta^-1 = (Theta J Theta^-1) U^-1 = J U^-1   if Theta and J commute.
# JU and J U^-1 are both involutions, and they commute iff U^4 = 1. So GENERICALLY the two
# projectors do NOT commute. Test that, then measure what survives.

set.seed(113)
n <- 6; N <- 2*n
Im_ <- diag(n); Z <- matrix(0,n,n)
J  <- rbind(cbind(Im_,Z), cbind(Z,-Im_))
sym <- function() { M <- matrix(rnorm(n*n),n,n); (M+t(M))/2 }
A  <- rbind(cbind(Z,sym()), cbind(sym(),Z))          # anticommutes with J, as in A.15
EV <- eigen(A); Uof <- function(t) Re(EV$vectors %*% diag(exp(EV$values*t)) %*% solve(EV$vectors))

cat("=== 1. measured facts only. My algebra slipped three times here.\n\n")
S <- rbind(cbind(Z,Im_), cbind(Im_,Z))
cat("   Attempt 1 assumed Theta and J commute and predicted Theta(JU)Theta^-1 = J U^-1.\n")
cat("   They anticommute here, so that was wrong. Attempt 2 corrected the sign to\n")
cat("   -J U^-1 and was ALSO wrong, because Theta M Theta^-1 = S M* S and S U S is not\n")
cat("   U^-1 unless S anticommutes with the generator, which it does not. A first pass\n")
cat("   also printed max|L2 - L2|, zero by construction, which checks nothing.\n")
cat("   Three slips on one chain. So: measurements only, and the decisive one in part 2\n")
cat("   needs none of this.\n\n")
cat(sprintf("   S and J anticommute:  |S J + J S| = %.2e\n\n", max(abs(S %*% J + J %*% S))))
cat("        t    Th(JU)Th^-1 an involution?   = JU?      = -JU^-1?    [JU,Th(JU)Th^-1]\n")
for (t in c(0.2, 0.7, 1.5)) {
  U <- Uof(t); L1 <- J %*% U
  C <- S %*% Conj(L1) %*% S
  cat(sprintf("   %8.2f %22.1e %12.2e %13.2e %18.2e\n", t,
      max(Mod(C %*% C - diag(N))), max(Mod(C - L1)), max(Mod(C + J %*% solve(U))),
      max(Mod(L1 %*% C - C %*% L1))))
}
cat("\n  Theta conjugates JU to a DIFFERENT involution, equal to neither JU nor -JU^-1,\n")
cat("  and the two do not commute. That is all part 2 needs.\n")

cat("\n=== 2. so how big is the common +1 eigenspace?\n\n")
cat("   Theta is ANTILINEAR, so its +1 eigenspace is a REAL subspace of real dimension N.\n")
cat("   JU is linear with a +1 eigenspace of real dimension N inside the same 2N-real\n")
cat("   dimensional space. Two N-dimensional subspaces of a 2N-dimensional space meet\n")
cat("   generically in dimension zero. Measure it.\n\n")
realify <- function(Mc) rbind(cbind(Re(Mc), -Im(Mc)), cbind(Im(Mc), Re(Mc)))
commondim <- function(t) {
  U <- Uof(t); L <- J %*% U
  # +1 eigenspace of the linear involution L, as a real subspace of C^N realified
  PL <- realify((diag(N)+0i + L)/2)
  # +1 eigenspace of Theta = K o S : vectors with S conj(v) = v. Realified projector:
  PT <- 0.5*rbind(cbind(diag(N) + S, matrix(0,N,N)), cbind(matrix(0,N,N), diag(N) - S))
  M <- rbind(PL - diag(2*N), PT - diag(2*N))
  sv <- svd(M)$d
  sum(sv < 1e-9) }
cat("        t      dim(+1 of JU)   dim(+1 of Theta)   dim of INTERSECTION\n")
for (t in c(0.2, 0.7, 1.5, 3.0)) 
  cat(sprintf("   %8.2f %15d %18d %21d\n", t, N, N, commondim(t)))
cat("\n  Zero. Generically no nonzero state satisfies both conditions.\n")

cat("\n=== 3. and when U^4 = 1, so the involutions DO commute?\n\n")
B <- rbind(cbind(Z, Im_*(pi/2)), cbind(-Im_*(pi/2), Z))   # generates a quarter turn: U^4 = 1
EVb <- eigen(B); Ub <- Re(EVb$vectors %*% diag(exp(EVb$values)) %*% solve(EVb$vectors))
cat(sprintf("   |U^4 - 1| = %.2e\n", max(abs(Re(EVb$vectors %*% diag(exp(4*EVb$values)) %*% solve(EVb$vectors)) - diag(N)))))
L1 <- J %*% Ub; L2 <- J %*% solve(Ub)
cat(sprintf("   |[JU, JU^-1]| = %.2e   -> they commute\n", max(abs(L1%*%L2 - L2%*%L1))))
cat("\n  So the compatible case is not empty, it is a condition on the dynamics: the\n")
cat("  evolution between the singular branches must satisfy U^4 = 1.\n")

cat("\n=== 4. flatly, and the caveat that may dissolve the whole question\n\n")
cat("  RESULT: the two projectors generically do not commute and their common fixed set\n")
cat("  is trivial. Imposed on the SAME data, A.15's singularity condition and the fold's\n")
cat("  Theta-invariance are incompatible unless U^4 = 1.\n\n")
cat("  CAVEAT, and it is serious: they may never apply to the same data. A.15's condition\n")
cat("  is on phase-space data at a BLACK HOLE singularity in the maximally extended\n")
cat("  solution; Theta-invariance is on the COSMOLOGICAL wavefunction. Those are different\n")
cat("  loci and different variables. The question only bites if the physical state is\n")
cat("  global - one state carrying both - and the paper does not settle whether it is.\n")
cat("  So this is a conditional incompatibility and it should be recorded as one.\n")
cat("  It is still worth having, because it is the first thing found that could force a\n")
cat("  choice BETWEEN two parts of the framework rather than adding to both.\n")

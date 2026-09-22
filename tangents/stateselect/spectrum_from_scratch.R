# Objection: the off-diagonal contributions cancel in N_- (giving
# N_- = diag(2P,1,1,2-2P) and Q = diag(P,1,1,2-P)), so the minimum eigenvalue is P, not n_*.
# Build the fermion operators explicitly with Jordan-Wigner and settle it. No algebra by hand.
# Basis order: |00>, |10>, |01>, |11>  (first slot = a mode, second = b mode)
I2 <- diag(2); Z <- diag(c(1,-1)); f <- matrix(c(0,1,0,0),2,2,byrow=TRUE)  # f|1>=|0>
a <- kronecker(f, I2)        # a acts on slot 1
b <- kronecker(Z, f)         # b acts on slot 2 with the JW string: THIS is the fermionic sign
ad <- t(a); bd <- t(b)

# check the algebra is actually fermionic before trusting anything built from it
ac <- function(X,Y) X%*%Y + Y%*%X
stopifnot(max(abs(ac(a,ad) - diag(4))) < 1e-12)
stopifnot(max(abs(ac(b,bd) - diag(4))) < 1e-12)
stopifnot(max(abs(ac(a,b)))  < 1e-12, max(abs(ac(a,bd))) < 1e-12)
cat("  anticommutators verified: {a,a+}={b,b+}=1, {a,b}={a,b+}=0\n")

Np <- ad%*%a + bd%*%b
cat(sprintf("  N_+ eigenvalues: %s  (trace %.1f)\n",
            paste(sprintf("%.3f", sort(Re(eigen(Np)$values))), collapse=" "), sum(diag(Np))))

for (P in c(0.9, 0.5, 0.2, 0.05, 0.01)) {
  s <- sqrt(P); cc <- sqrt(1-P)              # |s|^2 = P is the LZ transition probability
  am <- cc*a + s*bd
  bm <- cc*b - s*ad
  Nm <- t(am)%*%am + t(bm)%*%bm
  Q  <- (Np + Nm)/2
  offdiag <- max(abs(Nm - diag(diag(Nm))))
  ev <- sort(Re(eigen(Q)$values))
  nstar <- (1-sqrt(1-P))/2
  cat(sprintf("\n  P = %.2f\n", P))
  cat(sprintf("    largest off-diagonal element of N_-        : %.6f  %s\n", offdiag,
              ifelse(offdiag < 1e-12, "(claimed: they cancel)", "(they do NOT cancel)")))
  cat(sprintf("    eigenvalues of Q                           : %s\n",
              paste(sprintf("%.6f", ev), collapse="  ")))
  cat(sprintf("    trace of Q                                 : %.6f\n", sum(diag(Q))))
  cat(sprintf("    Claimed spectrum {P,1,1,2-P}          : %s\n",
              paste(sprintf("%.6f", sort(c(P,1,1,2-P))), collapse="  ")))
  cat(sprintf("    paper's claimed spectrum {n*,1/2,1/2,1-n*} : %s\n",
              paste(sprintf("%.6f", sort(c(nstar,.5,.5,1-nstar))), collapse="  ")))
  cat(sprintf("    min eigenvalue %.8f vs n_* = %.8f  -> %s\n", ev[1], nstar,
              ifelse(abs(ev[1]-nstar) < 1e-10, "MATCHES THE PAPER", "differs")))
}

cat("\n  ======================================================================\n")
cat("  RESOLUTION. Two separate things were wrong in the objection.\n")
cat("  ======================================================================\n")
cat("  (1) The off-diagonals of N_- do not cancel. With the Jordan-Wigner string on the\n")
cat("      b mode, the pair-creation terms ADD rather than cancel: a+b+ - b+a+ = 2a+b+,\n")
cat("      because the operators anticommute. Dropping the string makes them commute and\n")
cat("      the terms cancel, which is where diag(2P,1,1,2-2P) comes from. The table above\n")
cat("      shows off-diagonal elements of 0.20 to 1.00, so N_- is not diagonal.\n")
cat("  (2) 3.1 defines N_+- as the number PER PAIR MEMBER. The objection built the total\n")
cat("      number a+a + b+b, which has trace 4. Halving gives the paper's trace 2.\n\n")
cat("  Both corrections together reproduce the paper exactly:\n\n")
cat("        P        min eig(Q)/2        n_* = (1-sqrt(1-P))/2       difference\n")
for (P in c(0.95,0.9,0.8,0.7,0.64,0.5,0.2,0.05,0.01,0.001)) {
  s<-sqrt(P); cc<-sqrt(1-P)
  am <- cc*a + s*bd; bm <- cc*b - s*ad
  Q  <- (Np + t(am)%*%am + t(bm)%*%bm)/2
  half <- min(Re(eigen(Q)$values))/2
  ns <- (1-sqrt(1-P))/2
  cat(sprintf("   %8.3f %18.12f %26.12f %14.2e\n", P, half, ns, half-ns))
}
cat("\n  So the operator inequality Q >= n_* 1 stands as printed, and the ceiling of\n")
cat("  491.6 PeV that rests on it is unaffected.\n")

cat("\n  ======================================================================\n")
cat("  THE MATRIX, so that nothing has to be reconstructed.\n")
cat("  ======================================================================\n")
# NOTE THE CONVENTION. Np at the top of this file is the TOTAL number a+a + b+b, trace 4.
# The per-member operators 3.1 defines are half of that, so build them explicitly here rather
# than reusing Np. A first version of this block halved only the N_- term, which left
# Q[4,4] = 1.32 instead of 1 - P/2 and a determinant of 0.18 instead of P/4. The printed
# determinant check below is what caught it.
P <- 0.36; s <- sqrt(P); cc <- sqrt(1-P)
am <- cc*a + s*bd; bm <- cc*b - s*ad
Npm <- Np/2                                   # per pair member, out region
Nmm <- (t(am)%*%am + t(bm)%*%bm)/2            # per pair member, in region
Q  <- (Npm + Nmm)/2
cat(sprintf("  At P = %.2f, basis |00>,|10>,|01>,|11>:\n\n", P))
print(round(Q,6))
cat("\n  which is\n")
cat("      [ P/2                0    0    sqrt(P(1-P))/2 ]\n")
cat("      [ 0                 1/2   0         0         ]\n")
cat("      [ 0                  0   1/2        0         ]\n")
cat("      [ sqrt(P(1-P))/2     0    0     1 - P/2       ]\n")
cat(sprintf("\n  checked: Q[1,1]=%.8f vs P/2=%.8f;  Q[1,4]=%.8f vs sqrt(P(1-P))/2=%.8f\n",
            Q[1,1], P/2, Q[1,4], sqrt(P*(1-P))/2))
cat("\n  Pair creation mixes |00> with |11> and leaves the single-particle states alone, so\n")
cat("  the odd eigenvalues are exactly 1/2 and the even block has trace 1 and determinant\n")
cat(sprintf("  P/4 = %.8f (computed %.8f). Its eigenvalues are (1 +- sqrt(1-P))/2.\n",
            P/4, Q[1,1]*Q[4,4]-Q[1,4]^2))
cat("  That is a two-line derivation and it is now printed in 3.1, because three independent\n")
cat("  readings of the prose reconstructed three different wrong operators.\n")

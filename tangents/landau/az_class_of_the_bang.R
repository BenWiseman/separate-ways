# WHAT THIS FILE SHOWS, AND WHAT AN EARLIER VERSION WRONGLY CLAIMED.
#
# Earlier version: "the bang sits in Altland-Zirnbauer class CI, and the fold supplies the
# particle-hole half of that." Both halves are wrong, shown below by counterexample. What
# survives is a formal statement about a restricted matrix family and nothing about the
# physical cosmology.
I2 <- diag(2)
sx <- matrix(c(0,1,1,0),2)
sy <- matrix(c(0,1i,-1i,0),2)     # column-major: this IS [[0,-i],[i,0]]. Getting this backwards
sz <- diag(c(1,-1))               # cost an hour and nearly produced a false rebuttal.
anti <- function(U,M) U %*% Conj(M) %*% solve(U)
uni  <- function(U,M) U %*% M %*% solve(U)
sq   <- function(U) U %*% Conj(U)
z    <- function(M,tol=1e-12) max(Mod(M)) < tol
H    <- function(a,b,c=0,d=0) a*sz + b*sx + c*sy + d*I2

cat("=== WHAT HOLDS: the real traceless two-level family obeys the CI algebra ===\n")
a <- 0.7; b <- 0.35; Hp <- H(a,b)
cat(sprintf("  T = K     commutes  : %s    T^2 = %+d\n", z(anti(I2,Hp)-Hp), round(Re(sq(I2)[1,1]))))
cat(sprintf("  C = sy K  anticommut: %s    C^2 = %+d\n", z(anti(sy,Hp)+Hp), round(Re(sq(sy)[1,1]))))
cat(sprintf("  S = TC    anticommut: %s\n", z(uni(-sy,Hp)+Hp)))
cat("  (T^2,C^2) = (+1,-1) with chiral S. That is the CI algebra, for THIS matrix family.\n")

cat("\n=== WHAT FAILS 1: C^2 = -1 is supplied by TRACELESSNESS, not by the fold ===\n")
cat("  Theta covariance is Theta H(eta) Theta^-1 = H(-eta), i.e. sx conj(H(eta)) sx = H(-eta).\n")
cat("  Add an even identity term, which changes no physics the fold constrains:\n\n")
cat("     identity term   fold covariance defect   C anticommutation defect\n")
for (d in c(0, 0.2, 0.5)) {
  cat(sprintf("     %11.1f   %20.1e   %24.3f\n",
      d, max(Mod(anti(sx,H(a,b,0,d)) - H(-a,b,0,d))), max(Mod(anti(sy,H(a,b,0,d)) + H(a,b,0,d)))))
}
cat("\n  The fold is satisfied throughout while C fails. sy conj(H) sy = -H is automatic for any\n")
cat("  traceless Hermitian 2x2 matrix, so the particle-hole symmetry is a property of the\n")
cat("  truncation, not something the fold provides.\n")

cat("\n=== WHAT FAILS 2: an sy term does NOT drop the class to C ===\n")
cat("  T = K is lost, but a different antiunitary survives. With phi = arg(b + i c),\n")
cat("  T' = diag(e^-i phi, e^+i phi) K commutes with H at every point of the sweep:\n\n")
bb <- 0.35; cc <- 0.4; ph <- Arg(bb + 1i*cc)
Tp <- matrix(c(exp(-1i*ph),0,0,exp(1i*ph)),2,2)
Sp <- (cc*sx - bb*sy)/sqrt(bb^2+cc^2)
cat("        a        T' commutation defect\n")
for (aa in c(-2,-0.7,0,0.7,2))
  cat(sprintf("     %6.2f   %22.2e\n", aa, max(Mod(anti(Tp,H(aa,bb,cc)) - H(aa,bb,cc)))))
cat(sprintf("\n     T'^2 = +1 defect %.2e ;  S' = (c sx - b sy)/|.| anticommutes, defect %.2e\n",
    max(Mod(sq(Tp)-I2)), max(Mod(uni(Sp,H(0.7,bb,cc)) + H(0.7,bb,cc)))))
cat("  So the earlier claim that restoring sy costs time reversal and lands in class C is false.\n")

cat("\n=== WHAT FAILS 3: the cosmological map is not the T used here ===\n")
cat("  A.18's Theta sends a history at eta to one at -eta; T = K above is an INSTANTANEOUS\n")
cat("  operator at fixed eta. They are different objects and one does not imply the other.\n")
cat("  The many-body Theta^2 = (-1)^F likewise does not identify either single-particle AZ\n")
cat("  operator, and S(n) = S(1-n) is a permutation symmetry of a binary spectrum rather than\n")
cat("  an anticommutation of an antiunitary with a Hamiltonian.\n")

cat("\n  VERDICT. Retain only this: the real traceless two-level family admits operators obeying\n")
cat("  the CI symmetry algebra. We have NOT identified these with the physical symmetries of the\n")
cat("  cosmological field theory, and infer no AZ class and no random-matrix statistics for the\n")
cat("  bang. See NEXT_PAPERS.md F6, which is unaffected: RMT was already closed there for\n")
cat("  reasons independent of the class label.\n")

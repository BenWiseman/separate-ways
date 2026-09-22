# Named last turn: does the entanglement coefficient's symmetry under td -> 1/td survive a
# defect carrying a PHASE? That would separate "coincidence of the two-port structure" from
# "the inversion does real work".
#
# First: is the phase test even meaningful? On an OPEN chain a phase on a single bond is a
# pure gauge transformation, removable by redefining the operators beyond that bond. Check
# it rather than assume it, because if so the test I named is null and I need another.
#
# The test that DOES discriminate: c_eff is claimed to be a universal function of the
# TRANSMISSION. If so, a BOND defect and a SITE defect tuned to the same T must give the
# same coefficient, and the td -> 1/td symmetry is then inherited from T's symmetry rather
# than fundamental. If they differ, c_eff depends on more than T.

Sgen <- function(L, hmod) {
  h <- matrix(0+0i, L, L)
  for (i in 1:(L-1)) { h[i,i+1] <- -1; h[i+1,i] <- -1 }
  h <- hmod(h, L)
  e <- eigen(h); ord <- order(Re(e$values))[1:(L/2)]
  V <- e$vectors[, ord, drop=FALSE]
  C <- V %*% Conj(t(V))
  nu <- Re(eigen(C[1:(L/2),1:(L/2)])$values)
  nu <- pmin(pmax(nu,1e-14),1-1e-14); -sum(nu*log(nu)+(1-nu)*log(1-nu)) }
bond <- function(td, ph=0) function(h,L){ h[L/2,L/2+1] <- -td*exp(1i*ph); h[L/2+1,L/2] <- -td*exp(-1i*ph); h }
site <- function(eps)      function(h,L){ h[L/2,L/2] <- eps; h }
Tbond <- function(td) 4*td^2/(1+td^2)^2
Tsite <- function(eps) 4/(4+eps^2)
Ls <- c(64,128,256,512)
slope <- function(hmod) coef(lm(sapply(Ls, function(L) Sgen(L,hmod)) ~ log(Ls)))[2]

cat("=== 1. is a phase on the defect bond physical, or pure gauge?\n\n")
cat("        phase      S(L=200, td=0.7)       difference from phase 0\n")
s0 <- Sgen(200, bond(0.7, 0))
for (ph in c(0, pi/4, pi/2, pi)) 
  cat(sprintf("   %10.4f %20.6f %24.2e\n", ph, Sgen(200,bond(0.7,ph)), Sgen(200,bond(0.7,ph))-s0))
cat("\n  Identical to machine precision. On an open chain the phase IS pure gauge, so the\n")
cat("  test I named last turn cannot discriminate. Recorded rather than quietly dropped.\n")

cat("\n=== 2. the test that does: bond defect vs SITE defect at matched transmission\n\n")
s1 <- slope(bond(1.0))
cat("        target T     td (bond)    eps (site)    coeff/bond(T=1)   coeff/site(T=1)   agree?\n")
for (Tt in c(0.95, 0.88, 0.64, 0.30)) {
  td  <- uniroot(function(x) Tbond(x)-Tt, c(1e-4,1), tol=1e-12)$root
  eps <- uniroot(function(x) Tsite(x)-Tt, c(1e-6,50), tol=1e-12)$root
  a <- slope(bond(td))/s1; b <- slope(site(eps))/s1
  cat(sprintf("   %11.3f %12.5f %12.5f %17.5f %17.5f %8s\n", Tt, td, eps, a, b,
      if (abs(a-b) < 0.02) "yes" else "NO"))
}

cat("\n=== 3. reading it\n\n")
cat("  If the two columns agree, c_eff is a function of the transmission alone and the\n")
cat("  td -> 1/td symmetry is INHERITED from T's own symmetry under that map. The\n")
cat("  inversion would then be a property of how a two-port is parametrised, not a\n")
cat("  structure doing independent work, and the second appearance in A.13's rapidity\n")
cat("  form is the same fact twice rather than two facts.\n")
cat("  If they disagree, c_eff sees more than T and the question stays open.\n")

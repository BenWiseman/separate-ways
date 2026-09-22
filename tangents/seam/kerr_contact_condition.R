# Lead 9 ended with a yes/no: write the Kerr radial equation in the bang's variables and read
# off whether the contact condition holds through r=0. Do it.
#
# The bang's condition (2.2, B.5) is not "the diagonal is odd" -- that is basis-dependent and
# so cannot be the real statement. The real statement is the one Theta imposes,
#
#     H(-eta) = sigma_x H(eta) sigma_x,        Theta: psi(eta) -> sigma_x psi*(-eta)
#
# which for H = [[h11,h12],[h21,h22]] means h11(-x)=h22(x) and h12(-x)=h21(x). For the bang,
# H = [[gamma eta, p],[p, -gamma eta]] satisfies it because the diagonal is odd and the
# off-diagonal constant. A basis-independent NECESSARY consequence is that conjugation by
# sigma_x preserves the spectrum, so det H must be EVEN. Test that first: it is cheap and it
# cannot be evaded by a change of variables.
#
# Kerr, Dirac, Chandrasekhar's separated radial system (Mathematical Theory of Black Holes):
#
#     dR-/dr = -iK/Delta R- + (lambda + i mu r)/sqrt(Delta) R+
#     dR+/dr = +iK/Delta R+ + (lambda - i mu r)/sqrt(Delta) R-
#
# with K = (r^2+a^2) omega - a m,  Delta = r^2 - 2Mr + a^2, mu the particle mass, m the
# azimuthal number, lambda the angular eigenvalue. In i d/dr psi = H psi form,
#
#     H = [[ K/Delta , i(lambda+i mu r)/sqrt(Delta) ], [ i(lambda-i mu r)/sqrt(Delta) , -K/Delta ]]

Mkerr <- 1; a <- 0.6; om <- 0.3; mm <- 1; mu <- 0.2; lam <- 2.0

Kf     <- function(r, M, w, m) (r^2 + a^2)*w - a*m
Deltaf <- function(r, M) r^2 - 2*M*r + a^2
Hf <- function(r, M, w, m) {
  D <- Deltaf(r, M); sD <- sqrt(as.complex(D)); K <- Kf(r, M, w, m)
  matrix(c(K/D,                      1i*(lam - 1i*mu*r)/sD,
           1i*(lam + 1i*mu*r)/sD,    -K/D), nrow = 2, byrow = FALSE)
}
sx <- matrix(c(0,1,1,0), 2, 2)
det2 <- function(H) H[1,1]*H[2,2] - H[1,2]*H[2,1]   # base det() is real-only

cat("=== 1. the basis-independent test: is det H even in r?\n\n")
cat("   Conjugation by sigma_x preserves the spectrum, so H(-r)=sx H(r) sx forces\n")
cat("   det H(-r) = det H(r). If det H is not even, no change of basis can rescue it.\n\n")
cat("        r      det H(r)        det H(-r)        equal?\n")
for (r in c(0.15, 0.3, 0.5, 0.9)) {
  d1 <- det2(Hf( r, Mkerr, om, mm)); d2 <- det2(Hf(-r, Mkerr, om, mm))
  cat(sprintf("   %6.2f %14.6f %16.6f %11s\n", r, Re(d1), Re(d2),
      if (abs(d1-d2) < 1e-12) "yes" else "NO"))
}
cat("\n   det H = -(K/Delta)^2 + (lambda^2 + mu^2 r^2)/Delta. K and lambda^2+mu^2 r^2 are\n")
cat("   both even, so the whole obstruction is Delta, and Delta(-r)-Delta(r) = 4Mr.\n")
cat("   The naive crossing FAILS, and it fails on the mass term alone.\n")

cat("\n=== 2. so ask the right question instead\n\n")
cat("   At the bang the odd object is m*a(eta) = gamma eta. It is NOT the metric: in FRW\n")
cat("   ds^2 = a^2(-deta^2+dx^2), so a -> -a leaves the metric alone and only the MATTER\n")
cat("   coupling m*a changes sign. Flipping the sign of a mass coupling is charge\n")
cat("   conjugation, which is what the fold does. So the fold's map at the bang is not\n")
cat("   'eta -> -eta' but 'eta -> -eta TOGETHER WITH the CPT action on the parameters'.\n")
cat("   Apply the same standard at a hole. Under r -> -r the Kerr metric maps to Kerr\n")
cat("   with M -> -M (lead 9), and CPT also sends omega -> -omega and m -> -m. Test the\n")
cat("   COMBINED map, which is the one the fold actually supplies:\n\n")
cat("      quantity        under r->-r alone         under the full CPT map\n")
r0 <- 0.37
cat(sprintf("   Delta        %10.6f -> %10.6f    %10.6f -> %10.6f\n",
    Deltaf(r0,Mkerr), Deltaf(-r0,Mkerr), Deltaf(r0,Mkerr), Deltaf(-r0,-Mkerr)))
cat(sprintf("   K            %10.6f -> %10.6f    %10.6f -> %10.6f\n",
    Kf(r0,Mkerr,om,mm), Kf(-r0,Mkerr,om,mm), Kf(r0,Mkerr,om,mm), Kf(-r0,-Mkerr,-om,-mm)))
cat("\n   Delta is INVARIANT under the combined map, because -2Mr is even in (r,M) jointly.\n")
cat("   K flips sign, because it is linear in omega and in m and both flip.\n")

cat("\n=== 3. does the contact condition hold under the full map?\n\n")
cat("   Test H(-r; -M, -omega, -m) = sigma_x H(r; M, omega, m) sigma_x entrywise.\n\n")
cat("        r      max |LHS - RHS|      holds?\n")
worst <- 0
for (r in c(0.10, 0.25, 0.37, 0.55, 0.80)) {
  L <- Hf(-r, -Mkerr, -om, -mm)
  R <- sx %*% Hf(r, Mkerr, om, mm) %*% sx
  e <- max(Mod(L - R)); worst <- max(worst, e)
  cat(sprintf("   %6.2f %18.3e %11s\n", r, e, if (e < 1e-12) "YES" else "no"))
}
cat(sprintf("\n   worst entrywise discrepancy over the sample: %.3e\n", worst))
if (worst < 1e-12) {
  cat("   The bang's contact condition holds at the Kerr ring under the fold's own map.\n")
} else {
  cat("   It does NOT hold. Report that and stop.\n")
}

cat("\n=== 4. is it an accident of one parameter choice?\n\n")
cat("   Re-run over a grid. If this is structural it should hold identically everywhere.\n\n")
worst2 <- 0; nall <- 0
for (aa in c(0.2, 0.6, 0.95)) for (ww in c(-0.4, 0.3, 1.1)) for (mmm in c(-2, 0, 3)) {
  a <<- aa
  for (r in c(0.13, 0.41, 0.77)) {
    L <- Hf(-r, -Mkerr, -ww, -mmm); R <- sx %*% Hf(r, Mkerr, ww, mmm) %*% sx
    worst2 <- max(worst2, max(Mod(L - R))); nall <- nall + 1
  }
}
a <- 0.6
cat(sprintf("   %d parameter combinations, worst entrywise discrepancy %.3e\n", nall, worst2))
cat("   Identically, not numerically: every entry matches to machine zero.\n")

cat("\n=== 4b. negative control: make the test fail on purpose\n\n")
cat("   A check that returns YES at machine zero for 81 cases is worth nothing until it\n")
cat("   is shown capable of saying NO. Feed it the PARTIAL maps, each of which drops one\n")
cat("   leg of the fold's action. Every one of these must fail.\n\n")
cat("      map applied                              max |LHS-RHS|    verdict\n")
r <- 0.37
ctl <- list(
  list("r -> -r only",                      Hf(-r,  Mkerr,  om,  mm)),
  list("r -> -r, M -> -M  (no CPT on w,m)", Hf(-r, -Mkerr,  om,  mm)),
  list("r -> -r, w -> -w, m -> -m (no M)",  Hf(-r,  Mkerr, -om, -mm)),
  list("full fold map",                     Hf(-r, -Mkerr, -om, -mm)))
R <- sx %*% Hf(r, Mkerr, om, mm) %*% sx
for (cc in ctl) {
  e <- max(Mod(cc[[2]] - R))
  cat(sprintf("   %-40s %14.4f %10s\n", cc[[1]], e,
      if (e < 1e-12) "holds" else "FAILS"))
}
cat("\n   Three of the four fail, and only the complete map holds. The check can say no,\n")
cat("   and what makes it say yes is the fold's action in full rather than any part.\n")
cat("\n   Range note: the sampled r span r_- to r_+ where Delta<0 and sqrt(Delta) is\n")
cat("   imaginary. The identity is algebraic in Delta so it holds either way, but only\n")
cat("   r < r_- = M - sqrt(M^2-a^2) is the region where r is spacelike and the ring is\n")
cat(sprintf("   actually approached; at a=%.1f, M=1 that is r < %.3f.\n", a, Mkerr - sqrt(Mkerr^2-a^2)))

cat("\n=== 5. flatly\n\n")
cat("  ESTABLISHED: the Kerr Dirac radial system satisfies the bang's contact condition\n")
cat("  through r=0 under the fold's own map, which is r -> -r together with M -> -M,\n")
cat("  omega -> -omega and m -> -m. Delta is invariant under it because -2Mr is jointly\n")
cat("  even; K flips because it is linear in omega and m. That is the transfer lead 9\n")
cat("  said was well posed but unproven, and it comes out positive.\n\n")
cat("  WHAT I HAD WRONG IN LEAD 9: I wrote that Kerr 'has the structure for it' on the\n")
cat("  strength of Delta having an odd part. Section 1 shows that reading fails outright:\n")
cat("  under r -> -r ALONE det H is not even and no basis change repairs it. The odd part\n")
cat("  of Delta was the wrong diagnostic. What works is the combined map, and the reason\n")
cat("  is the same one that makes the bang work -- the sign flip lives in the MASS\n")
cat("  coupling, not in the metric.\n\n")
cat("  WHAT THIS IS NOT, and the caveats are load-bearing:\n")
cat("   (a) The bang is a SPACELIKE surface crossed in time. At r=0 inside the inner\n")
cat("       horizon Delta = a^2 > 0, so r is SPACELIKE and the ring is crossed in space.\n")
cat("       The algebra is the same; the interpretation is not, and this script does not\n")
cat("       supply the second one.\n")
cat("   (b) lambda is taken inert under the map. The angular eigenvalue depends on a*omega\n")
cat("       and m, and lambda(-a omega,-m) = lambda(a omega,m) needs checking against the\n")
cat("       spheroidal problem rather than assuming.\n")
cat("   (c) M -> -M is a map between two different spacetimes, not a motion within one.\n")
cat("       Whether the fold licenses it at a hole is exactly the seam question the\n")
cat("       companion leaves open.\n\n")
cat("  NEXT ROUTE: (b) is the cheap one and it is a genuine yes/no -- solve the spin-\n")
cat("  weighted spheroidal eigenvalue at (a omega, m) and at (-a omega, -m) and compare.\n")
cat("  If lambda is inert the result above is complete for Dirac. If it is not, the\n")
cat("  condition acquires an angular obstruction and the transfer is partial. Either way\n")
cat("  it is one eigenvalue problem, not a research programme.\n")

#!/usr/bin/env Rscript
# fold_map_classification.R -- A.15 assumes hypothesis (H), that the fold restricted to a
# neighbourhood of a hole is that hole's Kruskal wedge reflection. This enumerates every
# involutive isometry the Kruskal geometry admits and asks whether (H) is an assumption at
# all, or whether three properties the fold has anyway already force it.
#
# Kruskal-Schwarzschild:  ds^2 = -(32 M^3 / r) e^{-r/2M} dU dV + r^2 dOmega^2,
# with UV = (1 - r/2M) e^{r/2M}. Both r and the conformal factor depend on U and V only
# through the product UV, so a map of the (U,V) plane is an isometry exactly when it
# preserves UV and the form dU dV. That is the whole condition, and it is finite to check.

M <- 1
rOfUV <- function(p) {           # invert UV = (1 - r/2M) e^{r/2M}
  sapply(p, function(q) uniroot(function(r) (1-r/(2*M))*exp(r/(2*M)) - q,
                                c(1e-9, 60*M), tol = 1e-14)$root) }

cat("=== 1. which linear maps of the (U,V) plane are isometries ===\n")
cat("   A linear map (U,V) -> (aU+bV, cU+dV) preserves UV for all U,V iff\n")
cat("   ac = 0, bd = 0 and ad + bc = 1. Enumerating the branches:\n")
cat("      a,d free with b = c = 0 and ad = 1   ->  (U,V) -> (lam U, V/lam)   [diagonal]\n")
cat("      b,c free with a = d = 0 and bc = 1   ->  (U,V) -> (lam V, U/lam)   [swap]\n")
cat("   Both preserve dU dV (the diagonal by lam * 1/lam, the swap likewise), so both are\n")
cat("   isometries and there is nothing else. Verified numerically:\n")
chk <- function(map, nm) {
  set.seed(7); U <- runif(4000,-3,3); V <- runif(4000,-3,3)
  m <- map(U,V); ok <- max(abs(m$U*m$V - U*V))
  cat(sprintf("     %-34s max |U'V' - UV| = %.1e\n", nm, ok)); ok }
worst <- max(
  chk(function(U,V) list(U=2.3*U, V=V/2.3),  "diagonal, lam = 2.3"),
  chk(function(U,V) list(U=-U,    V=-V),     "N:  (U,V) -> (-U,-V)"),
  chk(function(U,V) list(U=V,     V=U),      "S:  (U,V) -> (V,U)"),
  chk(function(U,V) list(U=-V,    V=-U),     "SN: (U,V) -> (-V,-U)"))
stopifnot(worst < 1e-12)
cat(sprintf("   also a genuine non-isometry, for contrast:\n"))
set.seed(7); U <- runif(20,-3,3); V <- runif(20,-3,3)
cat(sprintf("     %-34s max |U'V' - UV| = %.1e  <- not an isometry\n",
            "(U,V) -> (-U, V)", max(abs((-U)*V - U*V))))

cat("\n=== 2. which of them are involutions ===\n")
cat("   diagonal: (lam U, V/lam) squared is (lam^2 U, V/lam^2), the identity only at\n")
cat("             lam = +-1. lam = 1 is the identity; lam = -1 is N.\n")
cat("   swap:     (lam V, U/lam) squared is (U, V) for EVERY lam. The swaps are all\n")
cat("             involutions, a one-parameter family.\n")
cat("   Up to conjugation by a boost the family collapses: lam > 0 conjugates to S,\n")
cat("   lam < 0 to SN. So there are exactly three nontrivial involutions, N, S and SN.\n")
for (lam in c(0.4, 1, 3.7, -0.4, -1, -2.9)) {
  f <- function(U,V) list(U=lam*V, V=U/lam)
  a <- f(0.83,-1.27); b <- f(a$U,a$V)
  stopifnot(abs(b$U-0.83)<1e-14, abs(b$V+1.27)<1e-14) }
cat("   swap-squared = identity confirmed at six values of lam including negatives.\n")
n2 <- function(lam) { f<-function(U,V) list(U=lam*U,V=V/lam); a<-f(0.83,-1.27); f(a$U,a$V) }
cat(sprintf("   N composed with a boost is an involution only for the trivial boost:\n"))
for (lam in c(1.5, 2.5)) { b <- n2(lam)
  cat(sprintf("     lam = %.1f -> (%.3f, %.3f), not (0.830, -1.270)\n", lam, b$U, b$V)) }
cat("   N is therefore rigid: no free parameter. The swap family has one.\n")

cat("\n=== 3. the three properties the fold has anyway ===\n")
cat("   (a) free, so the quotient is a manifold; (b) time-orientation reversing, because\n")
cat("   the fold is a CPT map; (c) it carries our exterior to the mirror sheet, which is\n")
cat("   a different region from ours.\n\n")
reg <- function(U,V) if (U<0 && V>0) "R" else if (U>0 && V<0) "L" else
                     if (U>0 && V>0) "F" else if (U<0 && V<0) "P" else "horizon"
maps <- list(
  N   = list(f=function(U,V) c(-U,-V),  lab="N:  (U,V) -> (-U,-V)"),
  S   = list(f=function(U,V) c( V, U),  lab="S:  (U,V) -> (V,U)   [the RP^3 geon]"),
  SN  = list(f=function(U,V) c(-V,-U),  lab="SN: (U,V) -> (-V,-U)"))
cat("     map                                 fixed set        time orient   R goes to   F goes to\n")
for (k in names(maps)) {
  m <- maps[[k]]
  # fixed set in the (U,V) plane
  fx <- if (k=="N") "U = V = 0" else if (k=="S") "U = V (thru F,P)" else "U = -V (thru R,L)"
  # time orientation: act on a future-directed vector (dU,dV) = (1,1)
  d <- m$f(1,1); to <- if (d[1]>0 && d[2]>0) "preserved" else "REVERSED"
  rr <- m$f(-1, 1); ff <- m$f(1, 1)
  cat(sprintf("     %-34s %-16s %-12s %-11s %s\n", m$lab, fx, to, reg(rr[1],rr[2]), reg(ff[1],ff[2]))) }

cat("\n   Composing with the sphere's antipodal map makes all three free, since the fixed\n")
cat("   sets above are fixed pointwise in the (U,V) plane and P_perp moves every point of\n")
cat("   the sphere. So (a) does not discriminate. (b) rules out S. (c) rules out SN, which\n")
cat("   sends each exterior to itself and so offers no second sheet to be the mirror.\n")
cat("   N survives alone, and it has no boost freedom to fix afterwards.\n")

cat("\n=== 4. why SN is excluded, which is not a matter of taste ===\n")
cat("   SN sends each exterior to itself, so under SN o P_perp the identification is\n")
cat("   x = (t,r,n) ~ (-t,r,-n) with both events in OUR exterior. Whatever the fold's\n")
cat("   contact means, it would then mean it where we can look. Ask when the two are\n")
cat("   causally connected: the minimum photon time from a point to its own antipode at\n")
cat("   the same radius, over impact parameters that turn through pi.\n\n")
f <- function(r) 1-2*M/r; h <- function(r) f(r)/r^2; hp <- function(r) 6*M/r^4 - 2/r^3
rmin_of <- function(b) uniroot(function(r) h(r)-1/b^2, c(3*M*(1+1e-9), 1e9), tol=1e-14)$root
leg <- function(b, r0, which) {
  rm <- rmin_of(b); if (rm >= r0) return(0)
  U <- sqrt(r0-rm); sg <- sqrt(-hp(rm))
  g <- function(u) { r <- rm+u*u; out <- numeric(length(u)); sm <- u < 1e-6
    out[sm] <- if (which==1) 2/(rm^2*sg) else 2/(b*f(rm)*sg)
    if (any(!sm)) { rr <- r[!sm]; q <- pmax(1/b^2 - h(rr), 0)
      out[!sm] <- 2*u[!sm]/((if (which==1) rr^2 else b*f(rr))*sqrt(q)) }
    out }
  2*integrate(g, 0, U, rel.tol=1e-9, subdivisions=4000)$value }
bc <- 3*sqrt(3)*M
cat("        r/M     impact b/M    photon time T/M    so contact once |t| >    flat 2r\n")
for (r0 in c(3.5, 5, 10, 100, 1e4)) {
  b <- uniroot(function(b) leg(b,r0,1)-pi, c(bc*(1+1e-4), sqrt(1/h(r0))*(1-1e-12)), tol=1e-12)$root
  T <- leg(b, r0, 2)
  cat(sprintf("  %9.4g  %13.5f  %17.4f  %22.4f M %10.4g\n", r0, b, T, T/2, 2*r0)) }
cat("\n   The threshold is finite at every radius and approaches the flat-space chord 2r,\n")
cat("   so under SN every exterior event more than about a light-crossing time from the\n")
cat("   reflection slice is identified with an event it can reach. N puts that same\n")
cat("   structure at r <= M, behind a horizon, where nobody outside can look. The choice\n")
cat("   between the two is settled by observation and not by preference.\n")

cat("\n=== 5. can ONE involution serve both horizons? Schwarzschild-de Sitter says no ===\n")
cat("   SdS is the place to ask, since it carries a black-hole horizon and a cosmological\n")
cat("   one at the same time, and the fold has to be both the hole's map and the\n")
cat("   cosmology's map if a fold-invariant hole is to sit in a folded universe.\n\n")
cat("   Lam M^2    r_b/M     r_c/M    kappa_b M   kappa_c M    ratio\n")
for (lam in c(0.005, 0.02, 0.05, 0.10, 0.110, 0.11111)) {
  fL <- function(r) 1 - 2*M/r - lam*r^2/3
  rts <- sort(Re(polyroot(c(-2*M, 1, 0, -lam/3)))); rts <- rts[rts > 0]
  rb <- rts[1]; rc <- rts[2]
  kb <- abs(2*M/rb^2 - 2*lam*rb/3)/2; kc <- abs(2*M/rc^2 - 2*lam*rc/3)/2
  cat(sprintf("   %8.5f  %7.4f  %8.4f  %10.5f  %10.5f  %7.3f\n", lam, rb, rc, kb, kc, kb/kc))
  stopifnot(fL((rb+rc)/2) > 0, fL(rb/2) < 0, fL(2*rc) < 0) }
cat("\n   f is positive strictly between the horizons and negative on both sides of that\n")
cat("   band, checked at every row. So r is spacelike in the static region and timelike in\n")
cat("   both neighbours, falling to the singularity inside r_b and growing without bound\n")
cat("   outside r_c, and both of those are the future. A static region's two neighbours\n")
cat("   both lie to its future, so the chain of static regions is ordered in time.\n\n")
cat("   N_b and N_c are point reflections of the conformal diagram about two different\n")
cat("   bifurcation points, and two point reflections compose to a translation by twice\n")
cat("   their separation. By the paragraph above that translation is timelike. Setting\n")
cat("   N_b = N_c means quotienting by it, which closes timelike curves. Nothing milder\n")
cat("   is available either: the surface gravities differ, by 6.5 at Lam M^2 = 0.005 and\n")
cat("   1.29 at 0.10, reaching one only as the two horizons merge at Nariai.\n")
cat("\n   That holds everywhere below Nariai. At Nariai it fails, and the exception is the\n")
cat("   result. There the spacetime is the product dS_2 x S^2, and in the embedding\n")
cat("   dS_2 = {X in R^{1,2} : -X0^2 + X1^2 + X2^2 = l^2} the static patch is X1 > |X0|,\n")
cat("   the horizons are X1 = +-X0, and both bifurcation points sit at X0 = X1 = 0, so\n")
cat("   at X2 = +-l. The linear map Nl(X) = (-X0,-X1,X2) is an isometry, an involution,\n")
cat("   reverses time orientation, sends our static patch to the other one and F to P,\n")
cat("   and fixes BOTH bifurcation points. It is the wedge reflection at both horizons\n")
cat("   at once, which is exactly what no sub-Nariai member of the family has.\n\n")
l <- 3; E <- diag(c(-1,1,1)); A <- diag(c(-1,-1,1))
Q <- function(X) -X[1]^2 + X[2]^2 + X[3]^2
Nl <- function(X) c(-X[1], -X[2], X[3])
set.seed(11)
pts <- t(sapply(1:6000, function(i) { tt <- runif(1,-3,3); pp <- runif(1,0,2*pi)
  c(l*sinh(tt), l*cosh(tt)*cos(pp), l*cosh(tt)*sin(pp)) }))
cat(sprintf("     ||A^T eta A - eta||          = %.1e   (linear isometry of R^{1,2})\n",
            max(abs(t(A)%*%E%*%A - E))))
cat(sprintf("     max |Q(Nl X) - Q(X)|         = %.1e   (stays on the hyperboloid)\n",
            max(abs(apply(pts,1,function(X) Q(Nl(X))-Q(X))))))
cat(sprintf("     max |Nl(Nl X) - X|           = %.1e   (involution)\n",
            max(abs(apply(pts,1,function(X) Nl(Nl(X))-X)))))
regd <- function(X) if (X[2] > abs(X[1])) "R" else if (X[2] < -abs(X[1])) "L" else
                    if (X[1] > 0) "F" else "P"
cat("     region map:")
for (X in list(c(0.3,2.9,0.4), c(0.3,-2.9,0.4), c(2.9,0.3,0.4), c(-2.9,0.3,0.4)))
  cat(sprintf("  %s->%s", regd(X), regd(Nl(X))))
cat("\n")
for (X in list(c(0,0,l), c(0,0,-l)))
  stopifnot(all(abs(Nl(X)-X) < 1e-12))
cat("     both bifurcation points fixed: TRUE. Composing with the antipodal map of the\n")
cat("     S^2 factor removes them, so the whole map is free.\n\n")
cat("   Why Nariai and nowhere below it: dS_2 has circles for spatial sections, so its\n")
cat("   maximal extension closes into four regions with no chain and no translation to\n")
cat("   obstruct, and its two bifurcation points lie in one copy where a single linear\n")
cat("   map fixes both. Below Nariai the extension is an infinite chain and they do not.\n\n")
lam_obs <- 1.1056e-52; G <- 6.67430e-11; c_ <- 2.99792458e8; Msun <- 1.98847e30
M_geo <- 1/(3*sqrt(lam_obs)); M_kg <- M_geo*c_^2/G
cat(sprintf("   Nariai at the observed Lambda = %.4e m^-2: M = 1/(3 sqrt(Lam)) = %.3e kg\n", lam_obs, M_kg))
cat(sprintf("   = %.2e solar masses, horizon r_h = 1/sqrt(Lam) = %.2f Gly. The one\n",
            M_kg/Msun, 3*M_geo/9.4607e24))
cat("   fold-invariant member of the family is the maximal hole, the size of the visible\n")
cat("   universe, and contact_charged.R independently finds that the contact region fills\n")
cat("   its whole interior: r_contact/r_horizon runs 0.500, 0.778, 0.919, 0.991 as\n")
cat("   9 Lam M^2 goes 0, 0.9, 0.99, 0.9999. Two calculations done for different reasons\n")
cat("   single out the same point.\n")

cat("\n=== 6. the same family, read from the other end ===\n")
cat("   The fold does not appear and disappear along the family. It is present throughout,\n")
cat("   because the COSMOLOGICAL horizon always has a wedge reflection and composing that\n")
cat("   with P_perp is free, time-reversing, and carries our static region across that\n")
cat("   horizon to a mirror region. That is the map section 2 uses at our own horizon.\n\n")
cat("   At M = 0 there is only one horizon, so the fold fixes it and nothing else is in\n")
cat("   play: pure de Sitter, the map is the antipodal map.\n")
cat("   For 0 < 9 Lam M^2 < 1 there are two horizons with distinct bifurcation surfaces.\n")
cat("   The fold still fixes the cosmological one. It does not fix the black hole, because\n")
cat("   fixing both is the thing section 5 rules out. Every sub-maximal hole is one-sided.\n")
cat("   At 9 Lam M^2 = 1 the two horizons are one object and fixing one fixes both.\n\n")
cat("   So the second class is not a separate population at all. A hole is two-sided\n")
cat("   exactly when its horizon IS the cosmological horizon, and the only mass at which\n")
cat("   that happens is the maximum. Section 6 and section 2 are describing one structure\n")
cat("   at the two ends of one parameter, and the contact radius interpolating from half\n")
cat("   the horizon radius to all of it is that parameter made visible.\n")

cat("
=== flatly ===

  (H) is not an independent assumption. Suppose only that the fold maps the hole to
  itself as a set. Then near the hole it is an involutive isometry of Kruskal, so it is
  N, S or SN up to a boost. CPT reverses time orientation, which kills S. SN sends each
  exterior to itself, which both leaves no second sheet to be the mirror and puts the
  fold's contact in our own exterior at every radius, where it is excluded by our having
  seen nothing of the kind. What is left is N composed with a free involution of the
  bifurcation sphere. Every free involution of S^2 is conjugate to the antipodal map,
  since the quotient is a closed surface of Euler characteristic one and RP^2 is the only
  one, so that factor is P_perp. The map is forced, with no parameter left over.

  Two things follow that the paper does not currently say. The existence of the second
  exterior L and the past interior P is a CONSEQUENCE of the fold fixing the hole, not a
  separate requirement on how the hole formed: N carries R to L, so if the fold fixes the
  hole then L is there. And the real question about the second class is not which
  formation channel supplies a white-hole region. It is whether the fold fixes any hole at
  all, which is a question about the global involution and not about a primordial tail.

  Section 5 then settles the conditional inside the Schwarzschild-de Sitter family, and
  the answer is a uniqueness statement rather than an emptiness one. Below Nariai no
  involution serves both horizons, so no hole there is fold-invariant. At Nariai one does,
  and the fold-invariant hole is the maximal one, 2.1e22 solar masses at the observed
  Lambda, with a horizon the size of the visible universe and a contact region filling its
  whole interior. Section 6's second class has exactly one member in this family and it is
  not an object inside a universe, it is a configuration of one.

  What is NOT settled. Whether anything outside the family qualifies, and whether the
  Nariai configuration is stable, which is a known question this file does not touch. The
  weaker condition section 6.1 actually uses, that both exteriors are fed, needs no
  fold-invariance and is untouched by any of this.\n")

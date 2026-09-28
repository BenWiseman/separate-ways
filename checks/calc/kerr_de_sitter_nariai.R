#!/usr/bin/env Rscript
# kerr_de_sitter_nariai.R -- does the fold-invariant class have rotating members?
#
# A.15 finds one fold-invariant hole in Schwarzschild-de Sitter, the Nariai one, and
# nariai_branch_charged.R widens that point to a curve with charge. Real black holes rotate,
# so the question that matters for Section 6 is whether spin does the same. This checks the
# two things that can be checked directly on the Kerr-de Sitter metric.
#
# Kerr-de Sitter, Boyer-Lindquist:
#   Delta_r = (r^2 + a^2)(1 - Lam r^2/3) - 2 M r,   Delta_th = 1 + (Lam a^2/3) cos^2(th),
#   Sigma = r^2 + a^2 cos^2(th),   Xi = 1 + Lam a^2/3.

lam <- 1/9
Dr   <- function(r,M,a,lm=lam) (r^2+a^2)*(1-lm*r^2/3) - 2*M*r
Dth  <- function(th,a,lm=lam) 1 + lm*a^2*cos(th)^2/3
Sig  <- function(r,th,a) r^2 + a^2*cos(th)^2
Xi   <- function(a,lm=lam) 1 + lm*a^2/3
gcomp <- function(r,th,M,a,lm=lam) {            # the five independent components
  S <- Sig(r,th,a); D <- Dr(r,M,a,lm); Dt <- Dth(th,a,lm); X <- Xi(a,lm); s2 <- sin(th)^2
  c(tt = -(D - Dt*a^2*s2)/(S*X^2),
    tp = -(Dt*(r^2+a^2) - D)*a*s2/(S*X^2),
    rr = S/D,
    hh = S/Dt,
    pp = (Dt*(r^2+a^2)^2 - D*a^2*s2)*s2/(S*X^2)) }

cat("=== 1. the transverse map is an exact isometry of Kerr-de Sitter ===\n")
cat("   Every metric function depends on theta through cos^2(theta) and sin^2(theta) and on\n")
cat("   phi not at all, so (theta, phi) -> (pi - theta, phi + pi) leaves all five alone. That\n")
cat("   is the same reason it works for Kerr, and Lambda does not disturb it.\n\n")
cat("      a        Lambda     worst |g(P x) - g(x)| over 3000 points\n")
set.seed(11)
for (p in list(c(0.3,1/9), c(0.8,1/9), c(1.2,1/9), c(0.8,0.02), c(0.8,0.30))) {
  a <- p[1]; lm <- p[2]
  rs <- runif(3000, 0.3, 12); ths <- runif(3000, 1e-3, pi-1e-3)
  e <- max(mapply(function(r,th) max(abs(gcomp(r,pi-th,1,a,lm) - gcomp(r,th,1,a,lm))), rs, ths))
  cat(sprintf("   %7.2f  %9.4f  %32.1e\n", a, lm, e))
  stopifnot(e < 1e-12) }
cat("\n   The check has to be able to fail, so plant a metric function odd in cos(theta):\n")
bad <- function(r,th,a) (r^2 + a^2*cos(th))
e <- max(sapply(runif(200,1,8), function(r) abs(bad(r,pi-1.1,0.8) - bad(r,1.1,0.8))))
cat(sprintf("     Sigma with cos instead of cos^2: worst difference %.3f   <- fails, as it must\n", e))
stopifnot(e > 1e-3)

cat("\n=== 2. the rotating Nariai locus, where the outer pair merges ===\n")
cat("   Eliminating M between Delta_r = 0 and Delta_r' = 0 gives the degenerate radius.\n")
cat("   Delta_r'' < 0 there means Delta_r <= 0 on both sides, so r is timelike and the\n")
cat("   near-horizon factor is dS_2, exactly as in the uncharged case.\n\n")
g <- function(r,a,lm=lam) r*(2*r*(1-lm*r^2/3) - 2*lm*r*(r^2+a^2)/3) - (r^2+a^2)*(1-lm*r^2/3)
d2 <- function(r,M,a,lm=lam,h=1e-5) (Dr(r+h,M,a,lm)-2*Dr(r,M,a,lm)+Dr(r-h,M,a,lm))/h^2
cat("      a       r0         M0       Delta_r''(r0)   inner horizon\n")
amax <- NA
for (a in c(0, 0.2, 0.5, 0.8, 1.0, 1.2, 1.28)) {
  r0 <- tryCatch(uniroot(function(r) g(r,a), c(2.0, 9), tol=1e-13)$root, error=function(e) NA)
  if (is.na(r0)) { cat(sprintf("   %7.2f   no outer merger\n", a)); next }
  M0 <- (r0^2+a^2)*(1-lam*r0^2/3)/(2*r0)
  z <- polyroot(c(a^2, -2*M0, 1-lam*a^2/3, 0, -lam/3))
  rr <- sort(Re(z[abs(Im(z))<1e-6])); inner <- rr[rr > -1e-9][1]
  cat(sprintf("   %7.2f %8.5f %10.5f  %+14.5f  %13.5f\n", a, r0, M0, d2(r0,M0,a), inner))
  stopifnot(abs(Dr(r0,M0,a)) < 1e-9, abs((Dr(r0+1e-6,M0,a)-Dr(r0-1e-6,M0,a))/2e-6) < 1e-6,
            d2(r0,M0,a) < 0)
  amax <- a }
cat(sprintf("\n   The branch runs from a = 0, where it is the Nariai point of the uncharged\n"))
cat(sprintf("   family, to a little past a = %.2f, beyond which the outer pair no longer\n", amax))
cat("   merges at this Lambda. Delta_r'' is negative along all of it, so the near-horizon\n")
cat("   factor is dS_2 and not AdS_2 at every spin, which is the condition that excluded the\n")
cat("   cold branch in the charged case.\n")

cat("\n=== 2b. where the rotating branch ends, in closed form ===\n")
cat("   The branch stops where the merged pair meets the inner horizon, so Delta_r has a\n")
cat("   TRIPLE root. The numbers say the quartic factorises there as\n")
cat("      Delta_r = -(Lam/3)(r - r0)^3 (r + 3 r0),\n")
cat("   and matching coefficients against (r^2+a^2)(1 - Lam r^2/3) - 2 M r gives\n")
cat("      Lam r0^2 = 2 sqrt(3) - 3,   a^2 = Lam r0^4,   M = (4 Lam/3) r0^3,\n")
cat("   so the spin at the endpoint is a pure number, free of Lambda:\n\n")
uu <- 2*sqrt(3) - 3
cat(sprintf("      a/M = 3/(4 sqrt(2 sqrt(3) - 3)) = %.10f\n\n", 3/(4*sqrt(uu))))
cat("      Lambda      r0 closed     r0 solved       a closed       a solved        a/M\n")
for (lm in c(1/9, 0.02, 0.5)) {
  gg  <- function(r,a) r*(2*r*(1-lm*r^2/3) - 2*lm*r*(r^2+a^2)/3) - (r^2+a^2)*(1-lm*r^2/3)
  top <- function(a) optimize(function(r) gg(r,a), c(0.2, 20/sqrt(lm*9)), maximum=TRUE, tol=1e-12)
  ac  <- uniroot(function(a) top(a)$objective, c(0.2/sqrt(lm*9), 8/sqrt(lm*9)), tol=1e-13)$root
  r0s <- top(ac)$maximum; M0 <- (r0s^2+ac^2)*(1-lm*r0s^2/3)/(2*r0s)
  r0c <- sqrt(uu/lm); acl <- sqrt(lm)*r0c^2
  cat(sprintf("   %10.5f  %11.6f  %11.6f  %13.6f  %13.6f  %10.6f\n",
              lm, r0c, r0s, acl, ac, ac/M0))
  stopifnot(abs(r0c-r0s) < 1e-6, abs(acl-ac) < 1e-6, abs(ac/M0 - 3/(4*sqrt(uu))) < 1e-6) }
cat("\n   Three values of Lambda, agreeing with the closed form to six figures. The rotating\n")
cat("   branch therefore runs from a = 0 to a/M = 1.10092 and stops, and that number is the\n")
cat("   same in every de Sitter background. Note it exceeds one: the endpoint is past the\n")
cat("   Kerr extremal bound, which a positive Lambda permits.\n")

cat("\n=== 3. what follows, and what does not ===\n")
cat("   Below the merger the obstruction of fold_map_classification.R section 5 applies\n")
cat("   unchanged: two horizons with distinct bifurcation surfaces, two wedge reflections,\n")
cat("   and a translation between them. At the merger there is one bifurcation surface and\n")
cat("   one reflection, so there is no translation to obstruct, and the transverse factor is\n")
cat("   handled by section 1 above at every spin. So the argument that gave one fold-invariant\n")
cat("   hole in the uncharged family gives a rotating one at every spin the branch reaches.\n\n")
cat("   NOT shown here. The rescaled near-horizon geometry of rotating Nariai is a warped\n")
cat("   product rather than a direct one, and this file does not write it down or repeat the\n")
cat("   embedding argument in it. What it establishes is that the two ingredients survive\n")
cat("   rotation: the transverse map is an exact isometry at every spin, and the radial\n")
cat("   factor is dS_2 along the whole branch. The uncharged endpoint is the case where the\n")
cat("   embedding argument was done in full.\n")

cat("
=== flatly ===

  The fold-invariant configuration is not isolated in spin either, and the branch has an
  exact end. The transverse map
  (theta, phi) -> (pi - theta, phi + pi) is an exact isometry of Kerr-de Sitter at every
  spin and every Lambda tested, for the same reason it works in Kerr: the metric functions
  depend on theta through cos^2(theta) alone. And the Nariai-type merger of the outer pair
  runs from a = 0 out to a little past a = 1.28 at Lambda = 1/9, with Delta_r'' negative
  along all of it, so the near-horizon radial factor is dS_2 and never AdS_2.

  The branch ends where Delta_r acquires a triple root, and there the quartic factorises as
  -(Lam/3)(r - r0)^3(r + 3 r0), which fixes Lam r0^2 = 2 sqrt(3) - 3 and gives

        a/M = 3/(4 sqrt(2 sqrt(3) - 3)) = 1.1009173688,

  the same in every de Sitter background and past the Kerr extremal bound, which a positive
  Lambda permits. Taken with the charged branch, Section 6's single member is better
  described as the boundary of the parameter space than as a point: maximal holes, where the
  black-hole and cosmological horizons coincide, carrying charge and spin up to that value. What is not done is the
  embedding argument inside the warped rotating geometry, which the uncharged case has in
  full and this one takes on structure rather than on calculation.\n")

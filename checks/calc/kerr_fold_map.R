#!/usr/bin/env Rscript
# kerr_fold_map.R -- REFEREE_OPEN item 6: which isometry the fold uses in Kerr. The answer
# is the same shape as in Schwarzschild, and the transverse factor is the part that has to
# be checked, because Kerr's bifurcation surface is not a round sphere.
#
# The fold is Theta = N o P_perp, with N the wedge reflection built from the horizon
# generator xi = d_t + Omega_H d_phi, and P_perp the antipodal map (theta, phi) ->
# (pi - theta, phi + pi). N lives in the Kruskal-type extension and fixes the bifurcation
# surface pointwise; everything that needs checking about P_perp can be checked in
# Boyer-Lindquist, and is checked here.

M <- 1
Sig <- function(r,th,a) r^2 + a^2*cos(th)^2
Del <- function(r,a)    r^2 - 2*M*r + a^2
g <- function(r,th,a) {                       # Boyer-Lindquist components
  S <- Sig(r,th,a); D <- Del(r,a); s2 <- sin(th)^2
  c(tt = -(1 - 2*M*r/S),
    tp = -2*M*a*r*s2/S,
    rr = S/D,
    hh = S,
    pp = (r^2 + a^2 + 2*M*a^2*r*s2/S)*s2) }

cat("=== 1. P_perp is an isometry of Kerr, at every spin ===\n")
cat("   Every metric function depends on theta only through cos^2(theta) and sin^2(theta),\n")
cat("   and none depends on phi, so (theta, phi) -> (pi - theta, phi + pi) leaves all five\n")
cat("   components alone. d(theta) changes sign and enters only squared.\n\n")
cat("      a/M      worst |g(P x) - g(x)| over 4000 points\n")
set.seed(5)
for (a in c(0.0, 0.3, 0.6, 0.9, 0.998)) {
  rp <- M + sqrt(M^2 - a^2)
  rs <- runif(4000, rp*1.001, 40); ths <- runif(4000, 1e-3, pi-1e-3)
  e <- max(abs(mapply(function(r,th) max(abs(g(r, pi-th, a) - g(r, th, a))), rs, ths)))
  cat(sprintf("   %7.3f   %28.1e\n", a, e)) }

cat("\n   The check must be able to fail, so try maps that are not isometries:\n")
bad <- list(`theta -> pi - theta only, in a metric with cos(theta)` =
              function(r,th,a) { S <- r^2 + a^2*cos(th); c(S, S) },
            `phi-dependence planted` = function(r,th,a) c(Sig(r,th,a)*(1+0.1*cos(th)), 0))
a <- 0.6
e1 <- max(sapply(runif(200,2,20), function(r) { th <- runif(1,0.2,pi-0.2)
        f <- bad[[1]]; max(abs(f(r,pi-th,a) - f(r,th,a))) }))
cat(sprintf("     a metric function linear in cos(theta): worst difference %.3f  <- fails\n", e1))
stopifnot(e1 > 1e-3)

cat("\n=== 2. P_perp is free, and it is the antipodal map ===\n")
cat("   Fixed points need theta = pi - theta, so theta = pi/2, and phi = phi + pi, which\n")
cat("   never holds. It is free on the whole spacetime, not only on the bifurcation\n")
cat("   surface, so Theta = N o P_perp is free wherever N is defined.\n")
cat("   Squaring: theta -> pi - (pi - theta) = theta and phi -> phi + 2pi = phi.\n")
cat("   On any 2-sphere of constant r and t it is the antipodal map, so the quotient of\n")
cat("   that sphere is RP^2, which is the same transverse structure the uncharged case has.\n")

cat("\n=== 3. what is different from Schwarzschild, and what is not ===\n")
cat("   Different: the bifurcation surface is not round. Its induced metric is\n")
cat("     ds^2 = Sigma_+ dtheta^2 + ((r_+^2 + a^2)^2 sin^2(theta)/Sigma_+) dphi^2,\n")
cat("   with Sigma_+ = r_+^2 + a^2 cos^2(theta), so the area radius varies with latitude.\n")
cat("   Not different: that metric still depends on theta through cos^2(theta) alone and\n")
cat("   not on phi, which is exactly what P_perp needs. Checking it separately:\n\n")
bif <- function(th, a) { rp <- M + sqrt(M^2-a^2); S <- rp^2 + a^2*cos(th)^2
                         c(hh = S, pp = (rp^2+a^2)^2*sin(th)^2/S) }
cat("      a/M     worst |h(P x) - h(x)| on the bifurcation surface    polar/equatorial    excess %\n")
for (a in c(0.0, 0.3, 0.6, 0.9, 0.998)) {
  ths <- runif(3000, 1e-3, pi-1e-3)
  e <- max(sapply(ths, function(th) max(abs(bif(pi-th,a) - bif(th,a)))))
  rp <- M + sqrt(M^2-a^2)
  ratio <- sqrt(bif(1e-6,a)["hh"]) / sqrt(bif(pi/2,a)["hh"])
  cat(sprintf("   %7.3f   %43.1e   %18.5f  %9.2f\n", a, e, ratio, 100*(ratio-1))) }
cat("\n   The surface gets less round as the spin rises, and P_perp does not care.\n")

cat("\n=== 4. the ambiguity A.15 flagged, and why it is not one ===\n")
cat("   A.15 declined to settle which isometry the fold uses, on the grounds that\n")
cat("   phi -> phi + pi and phi_K -> phi_K + pi are both isometries and differ by a time\n")
cat("   translation, so it was unclear which commutes with the wedge reflection. The\n")
cat("   ambiguity comes from defining the shift by a coordinate rather than by a flow.\n")
cat("   Define it as the flow of the axial Killing field psi = d_phi by parameter pi, and\n")
cat("   there is only one map and the commutation is immediate:\n\n")
cat("     N is built from the horizon generator xi = d_t + Omega_H d_phi and its Kruskal\n")
cat("     boost. Any isometry preserving xi preserves that structure and commutes with N.\n")
cat("       psi preserves xi, because [d_t, d_phi] = 0 and Omega_H is constant.\n")
cat("       E: theta -> pi - theta preserves d_t and d_phi separately, so it preserves xi.\n")
cat("     Hence Theta = N o exp(pi psi) o E satisfies Theta^2 = N^2 o exp(2 pi psi) o E^2 = 1,\n")
cat("     a genuine involution with nothing left to choose.\n\n")
cat("   Freeness does not even need N: exp(pi psi) o E fixes a point only if theta = pi/2\n")
cat("   and phi = phi + pi, and the second never holds. So Theta is free everywhere, and\n")
cat("   on the bifurcation surface, where N fixes every point, the freeness is carried\n")
cat("   entirely by the transverse factor, exactly as in the non-rotating case.\n\n")
cat("   Numerical corroboration of the two commutation claims, as metric invariances:\n")
a <- 0.7; rp <- M + sqrt(M^2-a^2); OmH <- a/(2*M*rp)
cat(sprintf("     a/M = %.1f, r_+ = %.5f, Omega_H M = %.6f\n", a, rp, OmH))
set.seed(9); rs <- runif(2000, rp*1.001, 30); ths <- runif(2000, 1e-3, pi-1e-3)
eE  <- max(mapply(function(r,th) max(abs(g(r,pi-th,a) - g(r,th,a))), rs, ths))
cat(sprintf("     E is an isometry:              worst |g(E x) - g(x)| = %.1e\n", eE))
cat("     exp(pi psi) is an isometry:    exact, no metric component depends on phi\n")
cat(sprintf("     xi is null on the horizon:     g_tt + 2 Om g_tp + Om^2 g_pp at r_+ = %.1e\n",
    max(sapply(runif(200,1e-3,pi-1e-3), function(th) { v <- g(rp,th,a)
      abs(v["tt"] + 2*OmH*v["tp"] + OmH^2*v["pp"]) }))))

cat("
=== flatly ===

  Item 6 is answered and the answer is the unsurprising one, which is worth knowing because
  the alternative was that rotation broke the construction. The fold in Kerr is the wedge
  reflection built from the horizon generator xi = d_t + Omega_H d_phi, composed with the
  transverse map (theta, phi) -> (pi - theta, phi + pi).

  The transverse map is an isometry of Kerr at every spin, exactly and not approximately,
  because the metric functions depend on theta only through cos^2(theta) and on phi not at
  all. It is free everywhere, it is an involution, and on the bifurcation surface it is the
  antipodal map, which is the same transverse structure the non-rotating case uses. The
  bifurcation surface is genuinely not round, its polar area radius exceeding its equatorial one
  by 1.17 per cent at a = 0.3M, 5.41 at 0.6M, 18.02 at 0.9M and 37.15 at 0.998M, and none of that
  matters to the map.

  It also removes the ambiguity A.15 recorded. The two candidate angle shifts differ only
  because one of them was named by a coordinate; defined as the flow of the axial Killing
  field, there is a single map, and it commutes with the wedge reflection because that field
  and the equatorial reflection both preserve the horizon generator. Theta squares to the
  identity with nothing left to choose.

  What this does NOT do: it does not carry the r <= M contact bound to Kerr. That bound came
  from an angular budget, and the budget on a distorted sphere is a separate calculation
  this file does not attempt. A.15's contact radius falling to about 0.89M at a = 0.9 stays
  a first pass and stays labelled one.\n")

# What the fold adds to the field equations, and why it cannot be absorbed into G.
#
# Running Clausius on the fold's own horizons puts the image stress into the balance alongside
# ordinary matter. If that stress were a local curvature polynomial it would renormalise Newton's
# constant and nothing else, and the construction would be general relativity with a shifted G.
# It is not, and the reason is worth the calculation: the image stress depends on the world
# function between a point and its FOLD IMAGE, which is a bi-scalar and not a local invariant.
#
# What that makes the construction is worth stating before the arithmetic. Not a modification of
# relativity with a new coupling, and not relativity either: relativity plus one term that carries
# no free parameter at all, since it is fixed by the geometry, the fold and the matter's departure
# from conformal invariance. And by the silence theorem that term vanishes outside every horizon.

cat("=== 1. the image stress on the one geometry where it is exact ===\n")
cat("   Einstein static universe, radius a, coupling xi. image_stress_conformal.R gives\n")
cat("        rho_img = (6 xi - 1) / ( 16 pi^2 a^2 (1 + cos eta) ),\n")
cat("   with eta = Delta t / a the separation to the image. Every component follows from it by\n")
cat("   tracelessness and homogeneity.\n")
a <- 1
rho_img <- function(eta, xi) (6 * xi - 1) / (16 * pi^2 * a^2 * (1 + cos(eta)))
cat("      eta      rho_img (xi = 0)   rho_img (xi = 1/6)   rho_img (xi = 1/3)\n")
for (e in c(0.5, 1.5, 2.5, 3.0)) {
  cat(sprintf("   %7.2f  %17.6f  %19.2e  %18.6f\n",
              e, rho_img(e, 0), rho_img(e, 1/6), rho_img(e, 1/3)))
}

cat("\n=== 2. the curvature of that geometry, which does not depend on eta at all ===\n")
cat("   R_tt = 0, R_ij = (2/a^2) g_ij, R = 6/a^2, and every higher invariant is a constant,\n")
cat("   because the Einstein static universe is homogeneous. So ANY local expression built from\n")
cat("   the metric and its derivatives is a constant in eta:\n")
cat(sprintf("      R_tt = %g      R = %g      R_ab R^ab = %g      Kretschmann = %g\n",
            0, 6 / a^2, 3 * (2 / a^2)^2, 12 / a^4))
cat("   (the spatial Ricci is 2/a^2 per direction and the sphere's Riemann gives 12/a^4)\n")

cat("\n=== 3. so no local form can fit it, and the fit says by how much it fails ===\n")
cat("   Try the most general local ansatz to second order, rho_local = alpha R_tt + beta R g_tt\n")
cat("   + gamma, which on this geometry is a CONSTANT whatever alpha, beta and gamma are. Least\n")
cat("   squares against rho_img over a window in eta:\n")
for (xi in c(0, 1/3)) {
  etas <- seq(0.4, 2.8, length.out = 40)
  y <- rho_img(etas, xi)
  best <- mean(y)                       # the best constant is the mean
  rms <- sqrt(mean((y - best)^2))
  cat(sprintf("      xi = %.3f:  best constant %.6f, residual RMS %.6f, which is %.0f%% of it\n",
              xi, best, rms, 100 * rms / abs(best)))
  stopifnot(rms / abs(best) > 0.3)
}
cat("   A local term cannot follow a function of the separation to the image, because the\n")
cat("   geometry does not know that separation. The image stress is NON-LOCAL.\n")

cat("\n=== 4. the plant: a genuinely local stress must be fitted perfectly by the same ansatz ===\n")
cat("   Hand the fitter a stress that IS local on this geometry, say alpha R + gamma:\n")
for (al in c(0.3, -1.7)) {
  y <- rep(al * 6 / a^2 + 0.2, 40)
  rms <- sqrt(mean((y - mean(y))^2))
  cat(sprintf("      alpha = %+.1f:  residual RMS %.2e, so the test passes what it should\n", al, rms))
  stopifnot(rms < 1e-12)
}

cat("\n=== 5. what that makes the field equations ===\n")
cat("   Clausius on the fold's horizons, with the image stress in the balance, returns\n")
cat("        G_ab + Lambda g_ab = 8 pi G ( T_ab^matter + T_ab^img[g, Theta] ),\n")
cat("   where the second source is a functional of the geometry and the fold rather than an\n")
cat("   independent field. Because it is not local it cannot be absorbed into G, so the\n")
cat("   construction is not general relativity with a shifted coupling. And it carries no free\n")
cat("   parameter: image_stress_conformal.R fixes its size by the matter's departure from\n")
cat("   conformal invariance and nothing else.\n")
cat("      conformally invariant matter:   T^img = 0 identically\n")
cat("      coupling xi away from 1/6:      T^img proportional to (1 - 6 xi)\n")
cat("      massive matter:                 the same slot in the Hadamard coefficient carries m^2\n")

cat("\n=== 6. where it is switched on, which is the part that decides the physics ===\n")
cat("   The silence theorem: the cross-sheet commutator vanishes identically outside any\n")
cat("   horizon, so T^img is zero everywhere an observation has ever been made. It is nonzero\n")
cat("   only inside a horizon, and there only inside the contact region, which is the inner\n")
cat("   half of the interior.\n")
M <- 1
cat(sprintf("      horizon at r = %.0f M, contact boundary at r = %.0f M, so the term lives on\n", 2*M, M))
cat("      the inner half of the interior and nowhere else.\n")
cat("   So the construction agrees with general relativity on every measurement ever made, by\n")
cat("   theorem rather than by tuning, and departs from it only where nothing has been measured.\n")
cat("   That is what makes it an alternative rather than a rival: there is no regime in which it\n")
cat("   has to explain away an agreement general relativity already has.\n")

cat("\n=== 7. the limits ===\n")
cat("   The exact image stress is known on the Einstein static universe and at a flat inversion,\n")
cat("   not on Schwarzschild, where the coefficient is still open. The non-locality argument does\n")
cat("   not depend on that: it needs only that T^img is a function of the separation to the image\n")
cat("   while local invariants are not, and the Einstein static universe exhibits that cleanly.\n")
cat("   Clausius applying on the fold's horizons as it does on Jacobson's was an assumption when\n")
cat("   this was written and is not one now: universality_local_fold.R derives it.\n")

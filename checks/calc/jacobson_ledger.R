# Can the fold supply Jacobson's two inputs? A ledger, and the kill of a route I opened
# the same day.
#
# Jacobson (1995) derives the Einstein equation from Clausius on local Rindler horizons and
# has to assume exactly two things: a temperature, and an entropy eta A proportional to area
# with eta universal. quarter_area_jacobson.R killed A.15's quarter as a source for eta and
# asked what the invariant half IS. contact_conjugacy_general.R answered that it is a
# CONJUGATE LOCUS and suggested that, since a conjugate locus is where the Van Vleck
# determinant diverges and therefore where an induced-gravity cutoff would have to sit, the
# fold draws its own cutoff. That suggestion is tested here before anything is built on it.
#
# 2026-10-06: section 4's verdict is the state of the argument when this was written. Two of its
# NOT-supplied lines have since been supplied: area_law_from_fold.R gives the entropy law from the
# entanglement of the state the fold requires, eta_is_universal.R gives the same eta at every
# horizon, and universality_local_fold.R gives horizons away from the fold's fixed locus. The
# section now says so after its verdict; nothing computed here changed.

M <- 1

cat("=== 1. how far the conjugate locus is from the horizon ===\n")
cat("   A UV cutoff has to be short. The contact surface sits at r = M and the horizon at\n")
cat("   r = 2M, and inside the horizon the radial direction is timelike, so the separation\n")
cat("   is a proper TIME: int_M^{2M} dr/sqrt(2M/r - 1).\n")
tau <- integrate(function(r) sqrt(r / (2 * M - r)), M, 2 * M, rel.tol = 1e-12)$value
cat(sprintf("   numerically %.8f M, closed form (pi/2 + 1) M = %.8f M\n", tau, pi/2 + 1))
stopifnot(abs(tau - (pi/2 + 1)) < 1e-7)
cat("   That is macroscopic. For a solar mass, in Planck units:\n")
lp_over_M_sun <- 1.616255e-35 / (1.4766e3)        # Planck length over the solar grav. radius GM/c^2
cat(sprintf("      (pi/2 + 1) M / l_P = %.3e\n", (pi/2 + 1) / lp_over_M_sun))
cat("   A cutoff 38 orders of magnitude above the Planck length does not renormalise G, and\n")
cat("   an entropy cut off there would be smaller than A/4G by the square of that, not equal\n")
cat("   to it. THE CONJUGATE LOCUS IS NOT A UV CUTOFF.\n")

cat("\n=== 2. and there is a second, worse problem: flat space has no conjugate locus ===\n")
cat("   Jacobson needs a local Rindler horizon at EVERY point in EVERY null direction, and\n")
cat("   locally every spacetime is flat. In flat space the Jacobi equation is J'' = 0, so\n")
cat("   J = lambda, the Van Vleck determinant is 1 everywhere, and there is no conjugate\n")
cat("   point on any geodesic at any separation. So at the horizons Jacobson actually uses\n")
cat("   the fold draws NO surface at all, let alone a short one.\n")
Jflat <- function(lam) lam
cat(sprintf("   Delta = lambda^2 / det J for a null pair in flat space, at lambda = 3: %.1f\n",
            3^2 / Jflat(3)^2))
cat("   The route is dead. Recorded as dead rather than left implied.\n")

cat("\n=== 3. what the fold does supply, checked ===\n")
cat("   (a) THE HORIZONS. The fold at a Rindler horizon is the wedge reflection composed\n")
cat("       with the transverse antipode, and the composition is frame-independent even\n")
cat("       though neither factor is. In any orthonormal frame\n")
cat("         J = diag(-1,-1,+1,+1) and P = diag(+1,+1,-1,-1),  J P = -Id,\n")
cat("       so a single map -Id is simultaneously the fold of every Rindler wedge through\n")
cat("       its fixed point, at every boost and every orientation. Checked on random frames:\n")
set.seed(11)
rand_frame <- function() {
  # a random Lorentz frame: boost by a random rapidity in a random direction, then rotate
  n <- rnorm(3); n <- n / sqrt(sum(n^2)); w <- runif(1, -1.5, 1.5)
  B <- diag(4); ch <- cosh(w); sh <- sinh(w)
  B[1, 1] <- ch; B[1, 2:4] <- -sh * n; B[2:4, 1] <- -sh * n
  B[2:4, 2:4] <- diag(3) + (ch - 1) * outer(n, n)
  Q <- qr.Q(qr(matrix(rnorm(9), 3, 3))); R <- diag(4); R[2:4, 2:4] <- Q
  B %*% R
}
worst <- 0
for (i in 1:200) {
  Fm <- rand_frame(); Fi <- solve(Fm)
  J <- Fm %*% diag(c(-1, -1, 1, 1)) %*% Fi         # wedge reflection in that frame
  P <- Fm %*% diag(c(1, 1, -1, -1)) %*% Fi         # transverse antipode in that frame
  worst <- max(worst, max(abs(J %*% P + diag(4))))
}
cat(sprintf("       worst |J P + Id| over 200 random boosts and rotations: %.2e\n", worst))
stopifnot(worst < 1e-10)
cat("       and the factors themselves are NOT frame-independent, or the check is empty:\n")
Fm <- rand_frame()
J1 <- Fm %*% diag(c(-1, -1, 1, 1)) %*% solve(Fm)
cat(sprintf("       one such J differs from diag(-1,-1,1,1) by %.3f\n",
            max(abs(J1 - diag(c(-1, -1, 1, 1))))))
stopifnot(max(abs(J1 - diag(c(-1, -1, 1, 1)))) > 0.1)

cat("\n   (b) THE TEMPERATURE. Already in the paper and not re-derived here: the cross-sheet\n")
cat("       correlator at a bifurcate Killing horizon is the direct one shifted by i beta/2,\n")
cat("       which fixes the two-sheet squeeze at tanh r = e^{-beta omega/2} with nothing\n")
cat("       left to choose. That is the thermofield double, and its temperature is Unruh's.\n")
bw <- c(0.5, 1, 2, 4)
cat("         beta omega    tanh r      implied T ratio\n")
for (b in bw) cat(sprintf("         %8.2f   %9.6f   %14.6f\n", b, exp(-b/2), 1))
cat("       (the ratio is 1 by construction; the content is that nothing is left free)\n")

cat("\n   (c) THE FOCUSING SURFACES. contact_conjugacy_general.R: contact is conjugacy, at\n")
cat("       every charge and in every dimension, and the caustic order is D-3. So the fold's\n")
cat("       own contact surfaces are null-focusing surfaces, which is the class of object\n")
cat("       Raychaudhuri and therefore Clausius operate on.\n")

cat("\n=== 4. the ledger, stated so the gap is visible ===\n")
cat("   supplied by the fold:   the horizons        (a), and frame-independently\n")
cat("   supplied by the fold:   the temperature     (b), with nothing left to choose\n")
cat("   supplied by the fold:   the focusing        (c), contact = conjugacy\n")
cat("   NOT supplied:           eta, the entropy per unit area\n")
cat("   NOT supplied:           universality at points off the fold's own fixed locus\n")
cat("\n   With eta assumed universal, Jacobson's argument then returns the Einstein equation\n")
cat("   with 8 pi G = 2 pi / eta. So the fold supplies the form and leaves the coupling: the\n")
cat("   field equations become a consequence of the fold plus one constant, rather than a\n")
cat("   consequence of the fold plus a temperature plus an entropy law plus a constant.\n")
cat("   That is the whole claim. It is one input short of a derivation and it is not\n")
cat("   presented as one.\n")
cat("   That was the ledger when this was written and is not the ledger now. The entropy law\n")
cat("   and its coefficient's universality are derived in area_law_from_fold.R and\n")
cat("   eta_is_universal.R, and the horizons off the fixed locus in universality_local_fold.R.\n")
cat("   What the fold does not supply is the value of eta, 1/4G, which is measured.\n")

cat("\n=== 5. what the fold's own stress does to that ledger ===\n")
cat("   image_stress_conformal.R found the fold's image stress proportional to (1 - 6 xi).\n")
cat("   Clausius balances the matter flux across the horizon against the area change, so an\n")
cat("   image stress enters on BOTH sides and its contribution is not automatically zero.\n")
cat("   For conformal matter the image term vanishes on the Einstein static universe but not\n")
cat("   at the flat inversion, where it is +1/(96 pi^2 R^4), so the statement available is that\n")
cat("   the fold's back-reaction on its own Clausius ledger is UNCOMPUTED, and it is the\n")
cat("   next quantity that matters rather than a detail. Its scale, at a horizon:\n")
for (R in c(1, 10, 100)) cat(sprintf("      R = %5.1f M   rho_img = %.3e / M^4\n", R, 1/(96*pi^2*R^4)))
cat("   falling as the fourth power, so it is a horizon-scale effect and not a global one.\n")

cat("\n=== 6. where the fold's own term starts to matter, which is not where expected ===\n")
cat("   Near contact the image null-null component is, from image_stress_conformal.R,\n")
cat("        T_kk = (1 - 6 xi)(1 - cos eta)/(8 pi^2 a^4 (1 + cos eta)^2)  ->  (1-6xi)/(pi^2 s^4)\n")
cat("   with s = a delta the PROPER distance off the contact surface, so it depends on the\n")
cat("   distance alone and not on the model's radius. Jacobson's balance is between\n")
cat("   8 pi G T_kk and R_kk, and the background R_kk is of order one over the curvature\n")
cat("   radius squared. Setting them equal:\n")
cat("        8 pi G (1-6xi)/(pi^2 s^4) ~ 1/a^2   =>   s ~ (8 (1-6xi)/pi)^{1/4} sqrt(l_P a).\n")
cat("   The fold's correction turns on at the GEOMETRIC MEAN of the Planck length and the\n")
cat("   curvature radius, not at the Planck length.\n\n")
lP <- 1.616255e-35                      # metres
rg <- function(Msun) 1.4766e3 * Msun     # GM/c^2 in metres
cat("      object                  M (solar)     a = GM/c^2        sqrt(l_P a)      in\n")
for (o in list(list("solar-mass hole", 1), list("Cygnus X-1", 21),
               list("GW150914 remnant", 62), list("M87*", 6.5e9), list("Sgr A*", 4.3e6))) {
  a <- rg(o[[2]]); s <- sqrt(lP * a)
  unit <- if (s < 1e-14) sprintf("%.2f fm", s * 1e15) else
          if (s < 1e-9)  sprintf("%.2f pm", s * 1e12) else sprintf("%.2e m", s)
  cat(sprintf("      %-22s %10.3g   %12.4g m   %12.3e m   %s\n", o[[1]], o[[2]], a, s, unit))
}
cat("\n   Nuclear for a stellar hole and atomic for a supermassive one, with the scale going\n")
cat("   as the square root of the mass. That is 19 orders of magnitude above the Planck\n")
cat("   length, so if the fold's correction to the field equations is real it is not hidden\n")
cat("   behind quantum gravity: it is hidden behind a horizon, which is a different kind of\n")
cat("   inaccessible and a weaker one.\n")
cat("   The scale check the other way: the correction at one Schwarzschild radius out is\n")
s1 <- sqrt(lP * rg(1)); a1 <- rg(1)
cat(sprintf("      s/a for a solar-mass hole = %.3e, so the surface is thin by that factor.\n", s1/a1))

cat("\n=== 7. the mass channel, which dominates and runs the other way ===\n")
cat("   The (1 - 6 xi) form is the massless non-conformally-coupled case, the one computed\n")
cat("   exactly. Real matter breaks conformal invariance with MASS, and in the Hadamard\n")
cat("   coefficient V_0 = Delta^{1/2}[m^2 + (xi - 1/6)R]/2 the mass sits in the same slot\n")
cat("   but carries different dimensions, giving T_kk ~ m^2/s^2 rather than 1/s^4. Redoing\n")
cat("   the balance 8 pi G m^2/s^2 ~ 1/a^2 with G = l_P^2 and m = 1/lambda_C:\n")
cat("        s ~ sqrt(8 pi) l_P m a = sqrt(8 pi) (l_P / lambda_C) a.\n")
cat("   That is LINEAR in the hole's radius and grows with the mass of the field, so the\n")
cat("   heaviest field sets the scale and heavier means the correction reaches FURTHER\n")
cat("   from the contact surface, not less far. The thickness as a fraction of the radius\n")
cat("   is sqrt(8 pi) l_P/lambda_C and does not depend on the hole at all.\n\n")
lam <- function(mev) 1.054571817e-34 / (mev * 1e6 * 1.602176634e-19 / 2.99792458e8)
cat("      field            m (MeV)   lambda_C (m)   s/a          s, solar hole   s, Sgr A*\n")
for (fld in list(list("electron", 0.511), list("muon", 105.66), list("tau", 1776.9),
                 list("bottom", 4180), list("top", 172570), list("Higgs", 125250))) {
  lc <- lam(fld[[2]]); frac <- sqrt(8 * pi) * lP / lc
  cat(sprintf("      %-12s %10.4g   %.4e   %.3e    %.3e m   %.3e m\n",
              fld[[1]], fld[[2]], lc, frac, frac * rg(1), frac * rg(4.3e6)))
}
cat("\n   The top quark sets it: about 0.1 pm inside a solar-mass hole and half a micron\n")
cat("   inside Sgr A*, against the xi channel's 0.15 fm. So the mass channel dominates by\n")
cat("   three orders and the scale is not Planckian by any reading.\n")
cat("   Flatly, the status: the 1/s^4 law is computed exactly, and the m^2/s^2 law is read\n")
cat("   off the slot m^2 occupies in V_0 rather than computed. The massive image sum at a\n")
cat("   caustic has not been done, and until it is, the table above is a scale and not a\n")
cat("   prediction. What is not in doubt either way is that the scale is nowhere near l_P:\n")
cat(sprintf("      ratio of the top-quark thickness to the Planck length, solar hole: %.2e\n",
            sqrt(8*pi) * lP / lam(172570) * rg(1) / lP))

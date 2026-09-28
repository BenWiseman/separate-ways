#!/usr/bin/env Rscript
# What the fold's own term would have to be worth to stop the focusing, and what stops it being
# a claim yet.
#
# Two results already sit in the release without being put beside each other. The image stress is
# identically zero outside every horizon and nonzero on the inner half of a black hole's interior
# in four dimensions and no others; and its null-null component there is negative, reached twice,
# once through the coupling slot of the Hadamard coefficient on an exactly solvable caustic and
# once through the mass slot at a Ricci-flat hole. A negative null-null stress is a failure of the
# null convergence condition, and the null convergence condition is the hypothesis Penrose's 1965
# theorem uses to turn a trapped surface into an incomplete geodesic. Nowhere in either manuscript
# is that said, so this file works out what can be said and what cannot.
#
# What can be said is arithmetic. For a null k the trace and Lambda terms drop out of the field
# equations, so R_ab k^a k^b = 8 pi G T_ab k^a k^b exactly and the sign of one is the sign of the
# other. A spherically symmetric ingoing radial null congruence in the Schwarzschild interior is
# shear-free, so Raychaudhuri reads theta' = -theta^2/2 - R_kk, and a negative R_kk works against
# the focusing rather than with it. Turning theta around needs R_kk < -theta^2/2, which at the edge
# of the contact region is a definite number: 8/r_h^2.
#
# What cannot be said is that the fold supplies it. The computed sign belongs to the contact
# direction, the null geodesic that turns through pi on the sphere, and the congruence Penrose's
# theorem follows is the radial one. Those are different null directions, and separating them needs
# the whole image stress tensor at a hole rather than one of its components. That tensor is the open
# item the magnitude fork is short of. So what follows is a threshold and a scaling, not a
# resolution, and the last section says exactly which number would close it.

TOL <- 1e-10
fail <- 0
note <- function(ok, what) {
  if (!ok) { cat("   *** FAILED:", what, "\n"); fail <<- fail + 1 }
  invisible(ok)
}
set.seed(4771)

cat("=== 1. for a null vector the field equations reduce to R_kk = 8 pi G T_kk ===\n")
cat("   G_ab + Lambda g_ab = 8 pi G T_ab contracted twice with a null k: g_ab k^a k^b = 0 kills\n")
cat("   both the Lambda term and the -R g_ab/2 in G_ab, so R_kk = 8 pi G T_kk with no residue.\n")
cat("   Checked on random Lorentzian metrics with random symmetric T and a k solved to be null:\n")
cat("      trial   |k|^2        R_kk - 8 pi G T_kk (from the field equations)\n")
for (tr in 1:4) {
  A <- matrix(rnorm(16), 4); g <- diag(c(-1, 1, 1, 1)) + 0.15 * (A + t(A)) / 2
  gi <- solve(g)
  Tm <- matrix(rnorm(16), 4); Tm <- (Tm + t(Tm)) / 2
  Lam <- rnorm(1); Gn <- 1 / (8 * pi)                       # so 8 pi G = 1
  # a null k: pick a spatial part, solve the quadratic for k^0
  ks <- rnorm(3)
  qa <- g[1, 1]; qb <- 2 * sum(g[1, 2:4] * ks); qc <- as.numeric(t(ks) %*% g[2:4, 2:4] %*% ks)
  k0 <- (-qb - sqrt(qb^2 - 4 * qa * qc)) / (2 * qa)
  k <- c(k0, ks)
  n2 <- as.numeric(t(k) %*% g %*% k)
  # read R_ab off the field equations: R_ab = 8 pi G T_ab + (Lambda - 8 pi G T/2 ... ) but the
  # contraction only needs that R_ab = 8 pi G T_ab + c g_ab for SOME scalar c, which is what
  # G_ab + Lambda g_ab = 8 pi G T_ab gives. Take c at random to show it cannot matter.
  cc <- rnorm(1)
  Rab <- 8 * pi * Gn * Tm + cc * g
  d <- as.numeric(t(k) %*% Rab %*% k) - 8 * pi * Gn * as.numeric(t(k) %*% Tm %*% k)
  cat(sprintf("   %7d %12.2e %42.2e\n", tr, n2, d))
  note(abs(n2) < 1e-10 && abs(d) < 1e-10, "the g_ab piece drops out on a null vector")
}
cat("   So nothing in the fold's extra term, in Lambda, or in the trace can change the sign: the\n")
cat("   null-null components of stress and of Ricci are locked together.\n")

cat("\n=== 2. the vacuum congruence focuses on its own, which is the thing to be stopped ===\n")
cat("   Ingoing radial null geodesics in the Schwarzschild interior, affine parameter with\n")
cat("   dr/dlambda = -1, have cross-sectional area 4 pi r^2 and so theta = -2/r. Spherical\n")
cat("   symmetry makes the shear vanish, and the vacuum makes R_kk vanish, so Raychaudhuri must\n")
cat("   read theta' = -theta^2/2 identically. Checked against the derivative of -2/r:\n")
cat("      r        theta      dtheta/dlambda     -theta^2/2      difference\n")
for (r in c(1.8, 1.0, 0.5, 0.2, 0.05)) {
  th <- -2 / r; dth <- -2 / r^2; rhs <- -th^2 / 2
  cat(sprintf("   %6.2f %10.4f %16.4f %15.4f %14.1e\n", r, th, dth, rhs, abs(dth - rhs)))
  note(abs(dth - rhs) < 1e-12, "the vacuum congruence satisfies Raychaudhuri exactly")
}
cat("   theta reaches -infinity at r = 0 with no help from any stress, which is the content of\n")
cat("   the singularity theorem in this one symmetric case.\n")

cat("\n=== 3. the threshold: how negative R_kk has to be to turn theta around ===\n")
cat("   theta' > 0 needs R_kk < -theta^2/2. On the contact region r <= r_h/2 the entry value is\n")
cat("   theta = -4/r_h, so the threshold at the edge is 8/r_h^2, and it grows as r^-2 inward:\n")
rh <- 2
cat("      r/r_h      theta        threshold |R_kk|\n")
for (f in c(0.5, 0.25, 0.1, 0.02)) {
  r <- f * rh; th <- -2 / r
  cat(sprintf("   %8.3f %11.3f %20.3f\n", f, th, th^2 / 2))
}
cat(sprintf("   At the edge of the contact region the threshold is %.4f and 8/r_h^2 is %.4f.\n",
            (2 / (rh / 2))^2 / 2, 8 / rh^2))
note(abs((2 / (rh / 2))^2 / 2 - 8 / rh^2) < TOL, "the edge threshold is 8/r_h^2")

cat("\n=== 4. two thresholds, one pointwise and one for a region of finite extent ===\n")
cat("   theta' = B - theta^2/2 with B = -R_kk >= 0 integrates in closed form, so the thresholds\n")
cat("   are algebra and the numerical work is a check on the algebra rather than the source of it.\n")
cat("   With a = sqrt(2B), and lambda measured from the entry at r_0 = r_h/2:\n")
cat("      B = 0            theta = theta_0/(1 + theta_0 lambda/2), reaching -infinity at\n")
cat("                       2/|theta_0|, which is the bound Penrose's argument uses\n")
cat("      |theta_0| > a    theta = a coth(a(lambda + C)/2), reaching -infinity at\n")
cat("                       (2/a) artanh(a/|theta_0|), later than the vacuum but still finite\n")
cat("      |theta_0| < a    theta = a tanh(a(lambda - lambda*)/2), never divergent, crossing\n")
cat("                       zero at lambda* = (2/a) artanh(|theta_0|/a)\n")
rh <- 2; r0 <- rh / 2; th0 <- -2 / r0; thr <- th0^2 / 2; L <- r0
lam_blow <- function(B) {
  if (B <= 0) return(2 / abs(th0))
  a <- sqrt(2 * B)
  if (a >= abs(th0)) return(Inf)
  (2 / a) * atanh(a / abs(th0))
}
theta_at <- function(B, lam) {
  if (B <= 0) return(th0 / (1 + th0 * lam / 2))
  a <- sqrt(2 * B)
  if (a < abs(th0)) { C <- (2 / a) * atanh(a / th0); return(a / tanh((a / 2) * (lam + C))) }
  lst <- (2 / a) * atanh(abs(th0) / a); a * tanh((a / 2) * (lam - lst))
}
integ <- function(B, s2 = 0, lam = 3, h = 1e-5, th_start = th0) {
  th <- th_start; l <- 0
  while (l < lam) {
    th <- th + h * (B - s2 - th^2 / 2); l <- l + h
    if (th < -1e7) return(list(lam = l, theta = th, conj = TRUE))
  }
  list(lam = l, theta = th, conj = FALSE)
}
cat(sprintf("   Entry: r_0 = %.3f, theta_0 = %.3f, theta_0^2/2 = %.3f, and 8/r_h^2 = %.3f\n",
            r0, th0, thr, 8 / rh^2))
note(abs(thr - 8 / rh^2) < TOL, "the pointwise threshold is 8/r_h^2")
cat("   (a) The vacuum saturates the bound, which is what makes it the marginal case:\n")
o <- integ(0)
cat(sprintf("       integrated %.5f, closed form %.5f, 2/|theta_0| = %.5f\n",
            o$lam, lam_blow(0), 2 / abs(th0)))
note(o$conj && abs(o$lam - lam_blow(0)) < 1e-2, "vacuum blow-up matches the closed form")
cat("   (b) Constant B over the whole run:\n")
cat("        B/thr    lambda_blow closed   integrated     theta(3) integrated   closed form\n")
for (m in c(0.5, 0.9, 0.99, 1.01, 2, 10)) {
  B <- m * thr; lb <- lam_blow(B)
  # run past the closed-form blow-up where there is one, so the comparison is not cut short
  o <- integ(B, lam = if (is.finite(lb)) min(1.3 * lb, 12) else 3)
  cat(sprintf("     %8.3f %20s %13s %21s %13s\n", m,
              if (is.finite(lb)) sprintf("%.5f", lb) else "none",
              if (o$conj) sprintf("%.5f", o$lam) else "none",
              if (o$conj) "-" else sprintf("%+.6f", o$theta),
              if (o$conj) "-" else sprintf("%+.6f", theta_at(B, 3))))
  if (m < 1) note(o$conj && abs(o$lam - lb) < 2e-2 && lb > 2 / abs(th0),
                  sprintf("B/thr = %.2f focuses later than vacuum, where the closed form says", m))
  # above threshold theta never diverges, and the integration must track the tanh solution,
  # which has NOT reached a by lambda = 3 near the threshold: comparing with a instead of with
  # the solution was the wrong check and it failed for the right reason.
  if (m > 1) note(!o$conj && abs(o$theta - theta_at(B, 3)) < 5e-4,
                  sprintf("B/thr = %.2f tracks the tanh solution", m))
}
lo <- 0.5 * thr; hi <- 3 * thr
for (it in 1:60) { mid <- (lo + hi) / 2; if (is.finite(lam_blow(mid))) lo <- mid else hi <- mid }
cat(sprintf("       bisection on the closed form: %.10f against 8/r_h^2 = %.10f\n",
            (lo + hi) / 2, 8 / rh^2))
note(abs((lo + hi) / 2 - 8 / rh^2) < 1e-9, "bisection finds 8/r_h^2")
cat("   (c) The region is finite, and that is the threshold that matters. Outside r < r_h/2 the\n")
cat("       image stress is zero, so the help stops at lambda = L. A theta still negative there\n")
cat("       reaches -infinity afterwards, at lambda = L + 2/|theta(L)|. So the congruence has to\n")
cat("       be brought to theta >= 0 INSIDE the region, and lambda* <= L is artanh(1/u) <= u with\n")
cat("       u = a L / 2:\n")
us <- uniroot(function(u) atanh(1 / u) - u, c(1.0001, 5), tol = 1e-13)$root
Bstar <- 2 * (us / L)^2
cat(sprintf("       u* = %.9f, B* = 2(u*/L)^2 = %.6f = %.4f/r_h^2 = %.4f x the pointwise value\n",
            us, Bstar, Bstar * rh^2, Bstar / thr))
note(abs(atanh(1 / us) - us) < 1e-10, "u* solves artanh(1/u) = u")
cat("        B/B*     theta(L)     zero reached inside?    otherwise -infinity at lambda =\n")
for (m in c(0.8, 0.99, 1.01, 1.2)) {
  B <- m * Bstar; a <- sqrt(2 * B); thL <- theta_at(B, L)
  lst <- (2 / a) * atanh(abs(th0) / a)
  cat(sprintf("     %8.3f %12.5f %22s %32s\n", m, thL, if (lst <= L) "yes" else "no",
              if (thL < 0) sprintf("%.5f", L + 2 / abs(thL)) else "never"))
  if (m < 1) note(thL < 0 && lst > L, sprintf("B = %.2f B* leaves the region converging", m))
  if (m > 1) note(thL > 0 && lst <= L, sprintf("B = %.2f B* leaves it expanding", m))
}
cat("       Confirmed by integrating inside the region and then with B switched off:\n")
for (m in c(0.8, 1.2)) {
  B <- m * Bstar; o1 <- integ(B, lam = L)
  o2 <- integ(0, lam = 6, th_start = o1$theta)
  cat(sprintf("        B = %.2f B*: theta(L) = %+.5f, then %s\n", m, o1$theta,
              if (o2$conj) sprintf("-infinity %.4f later", o2$lam) else "finite throughout"))
  note(if (m < 1) o2$conj else !o2$conj, sprintf("integration agrees at B = %.2f B*", m))
}
cat(sprintf("       So the number to beat is %.1f/r_h^2 and not 8/r_h^2, and the factor of %.3f\n",
            Bstar * rh^2, Bstar / thr))
cat("       between them is the price of the region being half the interior rather than all of it.\n")

cat("\n=== 5. what the fold supplies, and the number that is missing ===\n")
cat("   The mass slot's image stress diverges at the caustic as 1/D, where D is the distance to\n")
cat("   it: sign_through_the_mass_slot.R fits the power at D^-1.018 against the coupling slot's\n")
cat("   D^-2.005. A quantity that grows without bound crosses any finite threshold, so the\n")
cat("   question is not whether 11.5/r_h^2 can be reached but where. Writing |R_kk| = c m^2 / D with\n")
cat("   c a stress coefficient, the crossing of the finite-region threshold B* is at\n")
cat("   D* = c m^2 / B*. Be precise about which part of c is missing, because it is less than it\n")
cat("   looks. The Green function's amplitude at the caustic IS determined:\n")
cat("   caustic_amplitude_profile.R moves the Jacobi profile at fixed length and family volume on\n")
cat("   a non-symmetric surface of revolution and the heat kernel does not move, which leaves the\n")
cat("   affine length, c_1 = sqrt(pi) L, and at the Schwarzschild contact geodesic\n")
cat("   Delta^{1/2} -> 6.6092 M s^{-1/2}. What is not determined is the step from that Green\n")
cat("   function to the stress tensor's coefficient, which the fork ledger records as short by a\n")
cat("   Jacobian scaling. So the table below is parameterised in the part that is open and not in\n")
cat("   the part that is closed:\n")
cat("      c        m r_h        D* / r_h\n")
for (c0 in c(1e-2, 1, 1e2)) for (mr in c(1e-3, 1, 1e3)) {
  cat(sprintf("   %8.0e %10.0e %16.3e\n", c0, mr, c0 * mr^2 / (Bstar * rh^2)))
}
cat("   Those are ordinary numbers rather than Planckian ones, which is the point worth taking\n")
cat("   from the table: the threshold is not out of reach by orders of magnitude, it is a factor\n")
cat("   away in a quantity nobody has computed. Two things block the claim and both are named\n")
cat("   rather than buried. The computed sign belongs to the contact null direction and Penrose's\n")
cat("   congruence is the radial one, so the whole image stress tensor is needed and not one\n")
cat("   component; and a 1/D divergence at a caustic is exactly where a Hadamard expansion stops\n")
cat("   being trustworthy, so the scaling is a guide to where to look and not a value.\n")

cat("\n=== 6. plants ===\n")
cat("   (a) A POSITIVE R_kk must bring the conjugate point forward of the vacuum's lambda = r_0:\n")
o0 <- integ(0); op <- integ(-2 * thr)
cat(sprintf("       R_kk = 0: lambda = %.5f.   R_kk = +%.3f: lambda = %.5f\n",
            o0$lam, 2 * thr, op$lam))
note(op$conj && op$lam < o0$lam, "PLANT (a) fires: positive R_kk focuses sooner")
cat("   (b) The pointwise threshold written as |theta_0| rather than theta_0^2/2. At r_h = 2 the\n")
cat("       two coincide, so it is tested at r_h = 8 where |theta_0| = 0.5 and theta_0^2/2 = 0.125:\n")
rh8 <- 8; r08 <- rh8 / 2; t8 <- -2 / r08
i1 <- integ(0.13, th_start = t8, lam = 40 * r08)$conj
i2 <- integ(0.12, th_start = t8, lam = 40 * r08)$conj
cat(sprintf("       theta_0 = %.3f, theta_0^2/2 = %.4f, |theta_0| = %.3f\n", t8, t8^2/2, abs(t8)))
cat(sprintf("       B = 0.13 (past theta_0^2/2): conjugate point? %-5s\n", i1))
cat(sprintf("       B = 0.12 (short of it):      conjugate point? %-5s\n", i2))
note(!i1 && i2, "PLANT (b) fires: the threshold is theta_0^2/2 and not |theta_0|")
cat("   (c) Shear must raise the pointwise threshold by exactly sigma^2:\n")
a0 <- !integ(1.01 * thr, s2 = 0)$conj
a1 <- !integ(1.01 * thr, s2 = 1)$conj
a2 <- !integ(1.01 * thr + 1, s2 = 1)$conj
cat(sprintf("       sigma^2 = 0, B = 1.01 x pointwise: stopped = %s\n", a0))
cat(sprintf("       sigma^2 = 1, same B:               stopped = %s\n", a1))
cat(sprintf("       sigma^2 = 1, B raised by 1:        stopped = %s\n", a2))
note(a0 && !a1 && a2, "PLANT (c) fires: the threshold rises by exactly sigma^2")

cat("\n=== 7. superseded in part, and by which file ===\n")
cat("   image_stress_shell.R replaces the reading of the two thresholds above. They are correct\n")
cat("   about a CONSTANT stress and the stress is not constant: it diverges on the contact sphere,\n")
cat("   which is the OUTER edge of the region and therefore where an ingoing congruence enters, so\n")
cat("   the largest value arrives at entry and falls inward and the integrated effect at the edge\n")
cat("   is divergent. The numbers here survive as the LEVEL the real stress has to reach; what\n")
cat("   that file computes instead is the thickness of the shell in which it does.\n")

cat("\n=== 8. and superseded again, on the SIGN this time ===\n")
cat("   contact_and_radial_are_orthogonal.R changes what section 7's reading is worth. The fold joins\n")
cat("   a point to its image at fixed Schwarzschild t, so the contact geodesic is the E = 0 family\n")
cat("   and its tangent has no t component. A stress invariant under the interior\'s isometries gives\n")
cat("   T_kk = |f|(k^t)^2 (T^t_t - T^r_r) + r^2 (k^th)^2 (T^th_th - T^r_r), so the contact direction\n")
cat("   measures the anisotropy and the radial congruence measures T^t_t - T^r_r: orthogonal\n")
cat("   combinations, not two views of one number. Conservation with the trace gives\n")
cat("   T^t_t - T^r_r -> -2 (T^th_th - T^r_r) at the caustic, so the two carry opposite signs.\n")
cat("   So the negative value computed along the contact direction points to a POSITIVE radial one,\n")
cat("   and the reading below, that the fold removes Penrose\'s hypothesis, is not supported as it\n")
cat("   stands. One antecedent decides it and that file names it. Read section 9 as conditional.\n")

cat("\n=== 9. what to write, and what not to ===\n")
cat("   Write: the fold's extra term is nonzero only inside a horizon, only in four dimensions,\n")
cat("   and only on the inner half of the interior, and the sign computed there is the sign that\n")
cat("   the null convergence condition forbids. That condition is the hypothesis of the\n")
cat("   singularity theorems, so the one region where general relativity predicts its own failure\n")
cat("   is the one region where this construction differs from it, and the difference goes the\n")
cat("   way that removes the hypothesis. The threshold for it to matter is 8/r_h^2 pointwise at the\n")
cat("   edge of the contact region, and 11.5/r_h^2 once the region's finite extent is paid for.\n")
cat("   Do not write: that the singularity is resolved, or that the theorem is defeated. The\n")
cat("   component computed is not the component the theorem follows, and the tensor that would\n")
cat("   settle it is the open item.\n")
if (fail == 0) cat("\n   all checks passed\n") else {
  cat(sprintf("\n   %d CHECK(S) FAILED\n", fail)); quit(status = 1)
}

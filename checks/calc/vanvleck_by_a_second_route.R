#!/usr/bin/env Rscript
# The one typed number on fork 9's critical path, recomputed from the analytic tidal field.
#
# WHY. The calibration fork 9 has to hit is Delta^{1/2} -> 3.9004 M s^{-1/2}, and it is the
# product of an exact half and a computed one: c_1 = sqrt(pi) M(pi + 2) = 9.113236 M, exact, and
# the reduced Van Vleck factor Delta'^{1/2} = sqrt(lambda_tot / J2end) = 0.427995. J2end =
# 47.561945 is computed in contact_vanvleck.R from a Riemann tensor built NUMERICALLY, and then
# TYPED into eight other files: this one's consumers are caustic_amplitude_profile.R,
# amplitude_rule_meets_the_mode_sum.R, caustic_amplitude_form.R, caustic_length_invariance.R,
# caustic_replacement_rule.R, collective_jacobian_check.R, null_momentum_route.R and
# swept_length.R. Eight copies of a number that one file computes is eight chances to be stale,
# and it sits under both the calibration target and kappa, which sets the shell.
#
# So it is recomputed here from the closed-form tidal field, sharing no machinery with a
# numerically differentiated metric.
#
# THE GEODESIC. Take M = L = 1, horizon at r = 2. The contact geodesic is the E = 0 null family,
# so k^t = 0 and the curve lies in the (r, phi) surface at fixed t with
#   dr/dphi = +- r sqrt|f| = +- sqrt(r(2M - r)),   solved by   r = M(1 + sin phi),
# and phi from 0 to pi runs r from the contact sphere out to the horizon and back, which is
# Section 5's two legs. The affine parameter follows from dlambda = r^2 dphi / L, so
#   lambda_tot = int_0^pi (1 + sin phi)^2 dphi = pi + 4 + pi/2 = 3 pi/2 + 4,
# which is the affine length the release quotes, reached here in closed form.
#
# THE TIDAL FIELD. Schwarzschild is Ricci-flat, so the optical tidal matrix is trace-free and its
# two transverse eigenvalues are -+ 3 M L^2 / r^5: one direction focuses and the other spreads,
# which is the structure caustic_multiplicity.R relies on. In phi the Jacobi equation
# d^2 J/dlambda^2 +- (3 M L^2/r^5) J = 0 becomes d/dphi [J'/r^2] +- 3 J/r^3 = 0.

TOL <- 1e-9
fail <- 0
note <- function(ok, what) { if (!ok) { cat("   *** FAILED:", what, "\n"); fail <<- fail + 1 } }

M <- 1; L <- 1
rof <- function(p) M*(1 + sin(p))

cat("=== 1. the affine length, in closed form and by quadrature ===\n")
lam_closed <- 3*pi/2 + 4
n <- 2000000; ph <- pi*(seq_len(n) - 0.5)/n
lam_quad <- sum((1 + sin(ph))^2)*pi/n
cat(sprintf("   closed form 3 pi/2 + 4 = %.9f     quadrature = %.9f     diff %.1e\n",
            lam_closed, lam_quad, abs(lam_closed - lam_quad)))
note(abs(lam_closed - lam_quad) < 1e-8, "the affine length is 3 pi/2 + 4")
cat("   and the curve really is the E = 0 null geodesic: dr/dphi must equal sqrt(r(2M - r))\n")
d1 <- sapply(c(0.3, 0.9, 1.7, 2.6), function(p) M*cos(p))
d2 <- sapply(c(0.3, 0.9, 1.7, 2.6), function(p) { r <- rof(p); sign(cos(p))*sqrt(r*(2*M - r)) })
cat(sprintf("      max |dr/dphi - sqrt(r(2M-r))| over four angles: %.1e\n", max(abs(d1 - d2))))
note(max(abs(d1 - d2)) < 1e-12, "the contact curve solves the null condition")

cat("\n=== 2. the two Jacobi branches, integrated in phi ===\n")
# d/dphi [ J' / r^2 ] = -+ 3 J / r^3, written as a first-order pair with u = J'/r^2
run <- function(sgn, N = 4000000) {
  h <- pi/N; J <- 0; u <- 1                     # J(0) = 0, dJ/dlambda(0) = u(0) = 1
  f <- function(p, J, u) c(u*rof(p)^2, -sgn*3*M*L^2*J/rof(p)^3)
  p <- 0
  for (i in seq_len(N)) {
    k1 <- f(p, J, u)
    k2 <- f(p + h/2, J + h*k1[1]/2, u + h*k1[2]/2)
    k3 <- f(p + h/2, J + h*k2[1]/2, u + h*k2[2]/2)
    k4 <- f(p + h,   J + h*k3[1],   u + h*k3[2])
    J <- J + h*(k1[1] + 2*k2[1] + 2*k3[1] + k4[1])/6
    u <- u + h*(k1[2] + 2*k2[2] + 2*k3[2] + k4[2])/6
    p <- p + h
  }
  c(J = J, u = u)
}
foc <- run(+1); spr <- run(-1)
cat(sprintf("   focusing  branch: J(pi) = %14.9f   (must vanish: it is the caustic)\n", foc["J"]))
cat(sprintf("   spreading branch: J(pi) = %14.9f   against the release's 47.561945\n", spr["J"]))
note(abs(foc["J"]) < 1e-5, "the focusing branch vanishes at the far end, which is the caustic")
note(abs(spr["J"] - 47.561945) < 1e-4,
     "the spreading branch reproduces 47.561945 from the analytic tidal field")

cat("\n=== 3. so the reduced Van Vleck factor, and the calibration target ===\n")
Dv <- sqrt(lam_closed/spr["J"])
c1 <- sqrt(pi)*(pi + 2)*M
cat(sprintf("   Delta' = lambda_tot / J2end = %.6f,  Delta'^{1/2} = %.6f  (release 0.427995)\n",
            lam_closed/spr["J"], Dv))
cat(sprintf("   c_1 = sqrt(pi) M(pi + 2) = %.6f M\n", c1))
cat(sprintf("   Delta^{1/2} -> %.4f M s^{-1/2}   (release 3.9004)\n", Dv*c1))
note(abs(Dv - 0.427995) < 1e-5, "the reduced factor agrees to five decimals")
note(abs(Dv*c1 - 3.9004) < 1e-3, "and the calibration target reproduces")

cat("\n=== 4. the plants ===\n")
cat("   The tidal strength is 3 M L^2 / r^5 and nothing else. Scale it and the number must move:\n")
for (s in c(0.5, 0.9, 1.1, 2)) {
  runs <- function(sgn, N = 200000) {
    h <- pi/N; J <- 0; u <- 1; p <- 0
    f <- function(p, J, u) c(u*rof(p)^2, -sgn*s*3*M*L^2*J/rof(p)^3)
    for (i in seq_len(N)) {
      k1 <- f(p,J,u); k2 <- f(p+h/2,J+h*k1[1]/2,u+h*k1[2]/2)
      k3 <- f(p+h/2,J+h*k2[1]/2,u+h*k2[2]/2); k4 <- f(p+h,J+h*k3[1],u+h*k3[2])
      J <- J + h*(k1[1]+2*k2[1]+2*k3[1]+k4[1])/6
      u <- u + h*(k1[2]+2*k2[2]+2*k3[2]+k4[2])/6; p <- p + h
    }
    J
  }
  cat(sprintf("      tidal x %.1f:  spreading J(pi) = %10.4f,  focusing J(pi) = %10.4f\n",
              s, runs(-1), runs(+1)))
  note(abs(runs(-1) - 47.561945) > 1, sprintf("a tidal field scaled by %g does not give 47.56", s))
}
cat("   And the trace-free structure is what makes one branch vanish while the other spreads:\n")
cat("   they are the SAME equation with opposite sign, which is Ricci-flatness and not a choice.\n")

cat("\n=== 5. what this settles ===\n")
cat("   The one typed number under the calibration target is confirmed by a route sharing no\n")
cat("   machinery with the one that produced it: closed-form geodesic, closed-form affine length,\n")
cat("   analytic tidal field 3 M L^2 / r^5, against a numerically differentiated metric there.\n")
cat("   So fork 9 may use 3.9004 as its target. What is NOT settled is that eight files carry\n")
cat("   the number as a literal; a gate now checks they all agree with the file that computes it.\n")

if (fail == 0) cat("\n   all checks passed\n") else {
  cat(sprintf("\n   %d CHECK(S) FAILED\n", fail)); quit(status = 1)
}

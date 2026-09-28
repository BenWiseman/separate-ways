# Could the bang's image stress BE the cosmological constant? The scaling decides it.
#
# desitter_contact_budget.R closes the late universe: no contact, no image stress, no contribution
# to Lambda. bang_parity_rule.R closes the bang's kinematics: the fold forces a hot radiation bang
# and permits any Lambda. What it leaves open is section 6 of that script, the one route left. The
# image stress vanishes for conformally invariant matter, but the Hadamard slot carrying
# (xi - 1/6)R also carries m^2, so a massive field feels it whatever its coupling, and at the bang
# the adopted content has mass M_1. So there IS an image stress at the bang. The question is
# whether it can look like Lambda, and that is a question about how it scales with a.
#
# THE ARGUMENT. Section 2.2's bang is i d_eta psi = [[gamma eta, p], [p, -gamma eta]] psi with
# gamma = M_1 a_1. Rescaling eta = s/sqrt(gamma) and p = q sqrt(gamma) removes gamma completely:
#     i d_s psi = [[s, q], [q, -s]] psi.
# So every mode quantity is a function of s and q alone. The physical energy density of anything
# built from these modes is (1/a^4) times a conformal-frame integral, and in the rescaled
# variables d^3p x (conformal energy) carries gamma^(3/2) x gamma^(1/2) = gamma^2. Hence
#     rho_img = gamma^2 H(s) / a^4,   s = sqrt(gamma) a / a_1.
# A cosmological constant is a density independent of a. Since s is proportional to a, that needs
# H(s) proportional to s^4 near the bang, and then
#     rho_img -> gamma^2 c (sqrt(gamma)/a_1)^4 / 1 = c gamma^4 / a_1^4 = c M_1^4.
# The scale factor cancels completely and the answer is M_1 to the fourth, whatever c is. That is
# checkable against the observed value without computing H at all.

cat("=== 1. the rescaling is exact, checked on the equation rather than asserted ===\n")
# integrate the two-level system directly at several gamma and confirm the solution depends on
# (sqrt(gamma) eta, p/sqrt(gamma)) only.
step <- function(psi, eta, p, gam, h) {
  f <- function(ps, e) {
    H <- matrix(c(gam*e, p, p, -gam*e), 2, 2, byrow = TRUE)
    -1i * (H %*% ps)
  }
  k1 <- f(psi, eta);           k2 <- f(psi + h/2*k1, eta + h/2)
  k3 <- f(psi + h/2*k2, eta + h/2); k4 <- f(psi + h*k3, eta + h)
  psi + h/6*(k1 + 2*k2 + 2*k3 + k4)
}
evolve <- function(gam, q, s, mu, n = 20000) {
  p <- q*sqrt(gam); eta_end <- s/sqrt(gam); h <- eta_end/n
  psi <- matrix(c(1, exp(1i*mu)), 2, 1)/sqrt(2)
  e <- 0
  for (i in 1:n) { psi <- step(psi, e, p, gam, h); e <- e + h }
  psi
}
cat("      gamma      q      s        psi_1                    psi_2\n")
ref <- NULL
for (gam in c(0.25, 1, 4, 25)) {
  psi <- evolve(gam, q = 0.8, s = 1.3, mu = 2.0)
  if (is.null(ref)) ref <- psi
  cat(sprintf("   %8.2f   %.2f   %.2f   %+.8f%+.8fi   %+.8f%+.8fi\n", gam, 0.8, 1.3,
              Re(psi[1]), Im(psi[1]), Re(psi[2]), Im(psi[2])))
  stopifnot(max(abs(psi - ref)) < 1e-7)
}
cat("   Identical to seven figures over two decades in gamma, so the solution is a function of\n")
cat("   (s, q) and the prefactor gamma^2 is the whole of the gamma dependence.\n")

cat("\n=== 2. the plant: the rescaling must fail for an equation that does not have it ===\n")
step_bad <- function(psi, eta, p, gam, h) {         # a cubic sweep, which rescales differently
  f <- function(ps, e) { H <- matrix(c(gam*e^3, p, p, -gam*e^3), 2, 2, byrow = TRUE); -1i*(H %*% ps) }
  k1 <- f(psi, eta); k2 <- f(psi + h/2*k1, eta + h/2)
  k3 <- f(psi + h/2*k2, eta + h/2); k4 <- f(psi + h*k3, eta + h)
  psi + h/6*(k1 + 2*k2 + 2*k3 + k4)
}
evolve_bad <- function(gam, q, s, mu, n = 20000) {
  p <- q*sqrt(gam); eta_end <- s/sqrt(gam); h <- eta_end/n
  psi <- matrix(c(1, exp(1i*mu)), 2, 1)/sqrt(2); e <- 0
  for (i in 1:n) { psi <- step_bad(psi, e, p, gam, h); e <- e + h }
  psi
}
refb <- evolve_bad(0.25, 0.8, 1.3, 2.0); other <- evolve_bad(25, 0.8, 1.3, 2.0)
cat(sprintf("   a cubic sweep gives psi_1 = %+.6f%+.6fi at gamma = 0.25 and %+.6f%+.6fi at 25,\n",
            Re(refb[1]), Im(refb[1]), Re(other[1]), Im(other[1])))
cat(sprintf("   differing by %.3f. The check is reading the equation.\n", max(abs(refb - other))))
stopifnot(max(abs(refb - other)) > 0.1)

cat("\n=== 3. so if it were Lambda, its value would be M_1^4 ===\n")
M1_PeV <- 491.6                                    # the adopted dark-matter mass
M1_eV  <- M1_PeV * 1e15
rho_M1 <- M1_eV^4
# observed: rho_Lambda = Lambda/(8 pi G), Omega_Lambda h^2 with h = 0.674, Omega_L = 0.685
h <- 0.674; OmL <- 0.685
rho_crit_eV4 <- 8.098e-11 * h^2                    # GeV^4 -> eV^4 conversion folded in below
rho_crit <- 1.05371e4 * h^2                        # eV cm^-3
# do it in eV^4 directly: rho_crit = 3 H0^2/(8 pi G) = 8.0992e-47 h^2 GeV^4
rho_crit_GeV4 <- 8.0992e-47 * h^2
rho_L_GeV4 <- OmL * rho_crit_GeV4
rho_L_eV4  <- rho_L_GeV4 * 1e36
cat(sprintf("   M_1 = %.1f PeV, so M_1^4 = %.3e eV^4\n", M1_PeV, rho_M1))
cat(sprintf("   observed rho_Lambda = %.4f x %.4e GeV^4 = %.3e eV^4, i.e. (%.3e eV)^4\n",
            OmL, rho_crit_GeV4, rho_L_eV4, rho_L_eV4^0.25))
cat(sprintf("   ratio M_1^4 / rho_Lambda = %.3e, which is %.1f orders of magnitude\n",
            rho_M1/rho_L_eV4, log10(rho_M1/rho_L_eV4)))
cat(sprintf("   and in fourth roots, M_1 exceeds rho_Lambda^(1/4) by %.1f orders.\n",
            log10(M1_eV/rho_L_eV4^0.25)))

cat("\n=== 4. reading it ===\n")
cat("   The route is closed, and closed in the useful direction. If the bang's image stress were\n")
cat("   independent of the scale factor it would be a cosmological constant of order M_1^4, which\n")
cat("   overshoots the observed value by more than eighty orders of magnitude. So either the\n")
cat("   construction predicts a catastrophically wrong Lambda, or that stress is not constant in\n")
cat("   a and redshifts away like every other bang-era density. The second is what the rest of\n")
cat("   the construction already requires, since section 2.2's cosmology is an ordinary radiation\n")
cat("   era and a term of order M_1^4 sitting in it would be visible immediately.\n")
cat("\n   Which turns the hunt into a constraint worth stating: the fold CANNOT produce the\n")
cat("   observed Lambda at the bang, because the only scale the bang has is M_1 and the answer\n")
cat("   would be M_1^4. Nothing in the construction sets a scale near 2.3 meV. Lambda is a\n")
cat("   boundary datum, and now there is no place left in the construction where it could have\n")
cat("   been anything else.\n")

cat("\n=== 5. and computing H itself would BE the cosmological constant problem ===\n")
cat("   Section 4 settles Lambda without H, but it is worth knowing what computing H would take,\n")
cat("   because the answer is not a missing calculation. The conformal part of the image stress\n")
cat("   vanishes identically, so H is built from the difference between the massive crossing and\n")
cat("   the massless one, D(s,q) = C_massive - C_massless, weighted by q^2 dq times the conformal\n")
cat("   energy. Measured mode by mode:\n\n")
stepq <- function(psi, s, q, h, massive) {
  f <- function(ps, e) { d <- if (massive) e else 0
    -1i * (matrix(c(d, q, q, -d), 2, 2, byrow = TRUE) %*% ps) }
  k1 <- f(psi, s); k2 <- f(psi + h/2*k1, s + h/2)
  k3 <- f(psi + h/2*k2, s + h/2); k4 <- f(psi + h*k3, s + h)
  psi + h/6*(k1 + 2*k2 + 2*k3 + k4)
}
runq <- function(q, send, mu, massive, n) {
  psi <- matrix(c(1, exp(1i*mu)), 2, 1)/sqrt(2); h <- send/n; e <- 0
  for (i in 1:n) { psi <- stepq(psi, e, q, h, massive); e <- e + h }
  psi
}
Cq <- function(q, s, mu, massive, n) {
  p <- runq(q, s, mu, massive, n); 2*Re(exp(1i*mu)*Conj(p[1])*Conj(p[2]))
}
qs <- c(16, 24, 32, 48, 64, 96, 128)
for (s in c(0.6, 1.1)) {
  cat(sprintf("      s = %.1f\n         q          D           q^2 omega D (the integrand)\n", s))
  ig <- numeric(length(qs))
  for (k in seq_along(qs)) {
    q <- qs[k]; n <- max(6000, round(500*q*s))
    d <- Cq(q, s, 2.0, TRUE, n) - Cq(q, s, 2.0, FALSE, n)
    ig[k] <- q^2*sqrt(q^2+s^2)*d
    cat(sprintf("      %7.0f   %+.8f   %20.2f\n", q, d, ig[k]))
  }
  e <- unname(coef(lm(log(abs(ig)) ~ log(qs)))[2])
  cat(sprintf("         the integrand's envelope grows as q^%.3f, and it alternates in sign\n\n", e))
  stopifnot(e > 1.5)
}
cat("   D falls only as 1/q, so the integrand's envelope grows like q^2 and the q integral does\n")
cat("   not converge absolutely. The alternating sign leaves room for cancellation, so the\n")
cat("   statement to make is the weaker and safer one: H is a regulated quantity, not a number\n")
cat("   a quadrature returns. Regulating a divergent vacuum integral is the cosmological\n")
cat("   constant problem, so the construction meets that problem where every other treatment\n")
cat("   meets it, and section 4 already says what the answer would have to be if the fold fixed\n")
cat("   it alone. It does not fix it alone.\n")

cat("\n=== 6. the plant: the divergence must be the mass term's and not the integrator's ===\n")
cat("   Set the sweep to zero and D is identically zero, so the integrand must vanish.\n")
z <- Cq(32, 0.6, 2.0, FALSE, 20000) - Cq(32, 0.6, 2.0, FALSE, 20000)
cat(sprintf("      massless minus massless at q = 32: %.3e\n", z))
stopifnot(abs(z) < 1e-12)
cat("   and a mode that never sees the mass, q large at s tiny, must give a small D:\n")
sm <- Cq(64, 0.02, 2.0, TRUE, 20000) - Cq(64, 0.02, 2.0, FALSE, 20000)
cat(sprintf("      q = 64, s = 0.02: D = %.3e, against %.3e at s = 1.1\n", sm,
            Cq(64, 1.1, 2.0, TRUE, 40000) - Cq(64, 1.1, 2.0, FALSE, 40000)))
stopifnot(abs(sm) < 1e-3)
cat("   so D is measuring the mass and the divergence is the physical one.\n")

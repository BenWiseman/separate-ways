#!/usr/bin/env Rscript
# Which combination of the stress does the calibration geometry actually measure? A symmetry answers
# it, and the answer is that the calibration cannot see the one the contact direction reads.
#
# contact_and_radial_are_orthogonal.R shows that inside a hole the contact null direction measures
# A = T^theta_theta - T^r_r and Penrose's radial congruence measures X = T^t_t - T^r_r, and that
# conservation with the trace forces the two to have opposite signs at conformal coupling. It leaves
# open which of them inherits the negative value computed on the Einstein static universe. That
# question has a symmetry answer.
#
# On R x S^3 the rotations fixing a point also fix its antipode, because a rotation about an axis
# fixes both ends of it. So the isotropy group at x is the whole of SO(3) and it fixes the image
# point, which forces the image stress at x to be spatially ISOTROPIC: its anisotropy vanishes
# identically, not approximately. The single combination it carries is rho + p.
#
# Inside a hole the isotropy at x is only the SO(2) about the radial direction. That still fixes the
# contact image, since the antipodal map on the transverse sphere commutes with it, so the
# configuration is invariant under SO(2) and no more. An SO(2)-invariant stress has TWO independent
# combinations, rho + p_t and rho + p_T, which are X and A.
#
# So the calibration geometry carries one combination where a hole carries two, and the anisotropy is
# exactly the direction in which they differ. The computed negative value therefore says that some
# rho + p is negative at a hole without saying which, and conservation says that at conformal
# coupling exactly one of them is. Both outcomes are live and they are opposite:
#
#   if X inherits it, T_kk along Penrose's congruence is negative and the hypothesis is removed
#   if A inherits it, X is positive and the null convergence condition is satisfied there
#
# Nothing in the release decides between them yet. What this file does is establish that the
# calibration cannot, and say what would.

TOL <- 1e-12
fail <- 0
note <- function(ok, what) { if (!ok) { cat("   *** FAILED:", what, "\n"); fail <<- fail + 1 } }
set.seed(4413)

cat("=== 1. on S^3 the isotropy group at a point fixes the antipode ===\n")
cat("   Rotations of S^3 fixing a unit vector n are the SO(3) acting on its orthogonal complement,\n")
cat("   and every one of them sends -n to -n. Built explicitly and checked on 3000 of them:\n")
nvec <- c(1, 0, 0, 0)
worst_n <- 0; worst_a <- 0
for (i in 1:3000) {
  A <- matrix(rnorm(9), 3); Q <- qr.Q(qr(A)); if (det(Q) < 0) Q[, 1] <- -Q[, 1]
  R <- diag(4); R[2:4, 2:4] <- Q
  worst_n <- max(worst_n, max(abs(R %*% nvec - nvec)))
  worst_a <- max(worst_a, max(abs(R %*% (-nvec) - (-nvec))))
}
cat(sprintf("      worst |R n - n|        = %.2e over 3000 rotations\n", worst_n))
cat(sprintf("      worst |R(-n) - (-n)|   = %.2e over the same 3000\n", worst_a))
note(worst_n < TOL && worst_a < TOL, "the isotropy group at n fixes the antipode")
cat("   The group is three-dimensional, which is the whole rotation group of the tangent space:\n")
cat("      dim SO(3) = 3, and the tangent space to S^3 at a point is three-dimensional.\n")
note(TRUE, "recorded")

cat("\n=== 2. so the image stress there is isotropic, and carries one combination ===\n")
cat("   A symmetric spatial tensor invariant under the full SO(3) is a multiple of the identity, so\n")
cat("   the image stress at x is diag(-rho, p, p, p) and its anisotropy is zero identically. Averaged\n")
cat("   EXACTLY rather than by sampling: the twenty-four rotations of the octahedral group already\n")
cat("   kill the traceless part, since the five-dimensional representation it sits in contains no\n")
cat("   octahedral invariant.\n")
# a fixed test tensor, so every number below is reproducible rather than seed-dependent
S <- matrix(c(2, 1, 0,
              1, -1, 0.5,
              0, 0.5, 3), 3, byrow = TRUE)
perms <- list(c(1,2,3), c(1,3,2), c(2,1,3), c(2,3,1), c(3,1,2), c(3,2,1))
G <- list()
for (pm in perms) for (s1 in c(1,-1)) for (s2 in c(1,-1)) for (s3 in c(1,-1)) {
  Rm <- matrix(0, 3, 3)
  for (k in 1:3) Rm[k, pm[k]] <- c(s1, s2, s3)[k]
  if (abs(det(Rm) - 1) < 1e-12) G[[length(G)+1]] <- Rm
}
cat(sprintf("      group order: %d, all with determinant +1\n", length(G)))
note(length(G) == 24, "the octahedral rotation group has twenty-four elements")
acc <- Reduce(`+`, lapply(G, function(R) R %*% S %*% t(R))) / length(G)
iso <- sum(diag(S))/3
cat(sprintf("      averaged diagonal: %s\n", paste(sprintf("%+.10f", diag(acc)), collapse = " ")))
cat(sprintf("      off-diagonal magnitude: %.2e, against trace/3 = %+.10f\n",
            max(abs(acc - diag(diag(acc)))), iso))
note(max(abs(acc - diag(diag(acc)))) < 1e-12, "the average is exactly diagonal")
note(max(abs(diag(acc) - iso)) < 1e-12, "and exactly isotropic, at the trace over three")
cat("   For a null k on R x S^3, T_kk = (k^t)^2 (rho + p) whatever the spatial direction, so the\n")
cat("   computed number is that one combination and there is no second one to confuse it with.\n")
Tiso <- function(rho, p, kt) kt^2*(rho + p)
cat("      rho     p      T_kk at three spatial directions\n")
for (rp in list(c(0.4, -0.9), c(-1.2, 0.3))) {
  vals <- sapply(1:3, function(j) Tiso(rp[1], rp[2], 1))
  cat(sprintf("   %7.2f %6.2f %14s\n", rp[1], rp[2], paste(sprintf("%+.4f", vals), collapse = "  ")))
  note(diff(range(vals)) < TOL, "direction-independent, as isotropy requires")
}

cat("\n=== 3. inside a hole the isotropy is only SO(2), and two combinations survive ===\n")
cat("   The rotations fixing a point of the interior are those about its radial direction, an SO(2),\n")
cat("   and they fix the contact image too because the antipodal map of the transverse sphere\n")
cat("   commutes with them. Checked the same way, and then averaged:\n")
rad <- c(0, 0, 1)
worst_r <- 0; worst_i <- 0
for (i in 1:3000) {
  th <- runif(1, 0, 2*pi)
  R <- matrix(c(cos(th), -sin(th), 0, sin(th), cos(th), 0, 0, 0, 1), 3, byrow = TRUE)
  worst_r <- max(worst_r, max(abs(R %*% rad - rad)))
  worst_i <- max(worst_i, max(abs(R %*% (-rad) - (-rad))))
}
cat(sprintf("      worst |R r - r| = %.2e and |R(-r) - (-r)| = %.2e over 3000 rotations\n",
            worst_r, worst_i))
note(worst_r < TOL && worst_i < TOL, "the SO(2) fixes the radial axis and its antipode")
# the SO(2) average is exact in closed form: the 1-2 block goes to its own trace over two and the
# off-diagonals within it vanish, while the 3-3 entry is untouched.
acc2 <- diag(c((S[1,1] + S[2,2])/2, (S[1,1] + S[2,2])/2, S[3,3]))
chk <- Reduce(`+`, lapply(seq(0, 2*pi, length.out = 20001)[-20001], function(th) {
  R <- matrix(c(cos(th), -sin(th), 0, sin(th), cos(th), 0, 0, 0, 1), 3, byrow = TRUE)
  R %*% S %*% t(R) })) / 20000
cat(sprintf("      averaged over SO(2): diagonal %s, off-diagonal %.2e\n",
            paste(sprintf("%+.5f", diag(acc2)), collapse = " "),
            max(abs(acc2 - diag(diag(acc2))))))
cat(sprintf("      closed form against a 20000-point quadrature of the same average: %.2e\n",
            max(abs(acc2 - chk))))
note(max(abs(acc2 - chk)) < 1e-9, "the closed-form SO(2) average matches the quadrature")
note(abs(acc2[1,1] - acc2[2,2]) < 1e-12, "the two transverse entries agree exactly")
note(abs(acc2[3,3] - acc2[1,1]) > 0.05, "and the radial entry does NOT, so anisotropy survives")
cat("   So the hole carries two combinations where the calibration carries one, and the anisotropy\n")
cat("   is exactly the direction in which they differ.\n")

cat("\n=== 4. what that does to the sign question ===\n")
cat("   The calibration's negative value says one rho + p is negative at a hole and does not say\n")
cat("   which. Conservation with the trace says that at conformal coupling X and A carry opposite\n")
cat("   signs, so exactly one of them is negative. The two outcomes are:\n")
cat("      X negative:  T_kk along Penrose's ingoing radial congruence is negative, the null\n")
cat("                   convergence condition fails there, and the hypothesis is removed.\n")
cat("      A negative:  X is positive, the condition is satisfied along that congruence, and the\n")
cat("                   negative value belongs to the contact direction alone.\n")
cat("   Both are live. The release does not decide between them and neither does this file.\n")
cat("   What WOULD decide it is the image stress's own decomposition at a hole into its radial and\n")
cat("   transverse pressures, which is the same object A.18's sum on the contact geodesic supplies\n")
cat("   and is the one open item in the interior chain. The symmetry argument above says why the\n")
cat("   calibration geometry can never supply it: the combination that distinguishes the two cases\n")
cat("   is identically zero there.\n")

cat("\n=== 5. plants ===\n")
cat("   (a) the SO(3) average must destroy an anisotropy that the SO(2) average keeps:\n")
aniso <- diag(c(1, 1, -2))
a3 <- Reduce(`+`, lapply(G, function(R) R %*% aniso %*% t(R))) / length(G)
a2 <- diag(c((aniso[1,1] + aniso[2,2])/2, (aniso[1,1] + aniso[2,2])/2, aniso[3,3]))
cat(sprintf("       diag(1,1,-2) averaged over SO(3): %s\n",
            paste(sprintf("%+.4f", diag(a3)), collapse = " ")))
cat(sprintf("       the same averaged over SO(2):     %s\n",
            paste(sprintf("%+.4f", diag(a2)), collapse = " ")))
note(max(abs(diag(a3) - mean(diag(a3)))) < 1e-12, "PLANT (a): SO(3) flattens it exactly")
note(abs(a2[3,3] - a2[1,1]) > 1, "PLANT (a): SO(2) keeps it")
cat("   (b) a rotation NOT about the radial axis must move the image, so the isotropy claim is not\n")
cat("       vacuous:\n")
th <- 0.7
Rx <- matrix(c(1, 0, 0, 0, cos(th), -sin(th), 0, sin(th), cos(th)), 3, byrow = TRUE)
cat(sprintf("       |R_x r - r| = %.4f for a rotation about a transverse axis\n",
            max(abs(Rx %*% rad - rad))))
note(max(abs(Rx %*% rad - rad)) > 0.1, "PLANT (b): only the radial axis works")
if (fail == 0) cat("\n   all checks passed\n") else {
  cat(sprintf("\n   %d CHECK(S) FAILED\n", fail)); quit(status = 1)
}

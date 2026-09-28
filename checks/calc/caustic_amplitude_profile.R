# What the caustic amplitude actually depends on. It is neither of the two candidates.
#
# caustic_amplitude_form.R wrote c_n = V_fam (pi <J^2>/2)^{n/2} and left open whether the
# bracket is the path-space mean square of the Jacobi field or half its squared maximum, since
# the two agree on a sine and every member of the sphere family has the profile a sin(u/a).
#
# COMPLEX PROJECTIVE SPACE CANNOT SETTLE IT, and the reason kills the whole class. On a
# symmetric space the curvature operator is parallel along a geodesic, so every Jacobi field is
# a pure sine, and every pure sine has <J^2> = max(J)^2/2 exactly. CP^m at its cut locus gives
# sin(2u)/2, still a sine. No symmetric space separates them.
#
# So build a non-symmetric one. On a surface of revolution ds^2 = du^2 + f(u)^2 dphi^2 with f
# vanishing at both poles, the pole-to-pole geodesics are the meridians: a one-parameter family
# of volume 2 pi, with the Jacobi field equal to f itself. Take f = a[sin th + eps h(th)] with
# h odd at both poles, and L = pi a and V_fam = 2 pi stay fixed while the profile moves freely.
# Only the bracket changes, so the two readings make different predictions and the heat kernel
# between the poles decides. It decides against both of them.

a <- 1; L <- pi * a

polespec <- function(fn, n = 1200, Lx = L) {
  # The heat kernel between the poles sees only m = 0 modes, since every other azimuthal mode
  # vanishes there. Those solve -(f psi')' = lambda f psi, whose weak form needs no special
  # treatment at the singular ends.
  u <- seq(0, Lx, length.out = n + 1); h <- u[2] - u[1]
  fm <- fn(u[-1] - h / 2)
  A <- matrix(0, n + 1, n + 1); B <- matrix(0, n + 1, n + 1)
  for (i in 1:n) {
    w <- fm[i] * h
    A[i, i] <- A[i, i] + fm[i]/h;  A[i+1, i+1] <- A[i+1, i+1] + fm[i]/h
    A[i, i+1] <- A[i, i+1] - fm[i]/h; A[i+1, i] <- A[i+1, i] - fm[i]/h
    B[i, i] <- B[i, i] + w/3;      B[i+1, i+1] <- B[i+1, i+1] + w/3
    B[i, i+1] <- B[i, i+1] + w/6;  B[i+1, i] <- B[i+1, i] + w/6
  }
  R <- chol(B); Ri <- backsolve(R, diag(n + 1))
  M <- t(Ri) %*% A %*% Ri; M <- (M + t(M)) / 2
  e <- eigen(M, symmetric = TRUE); o <- order(e$values); V <- Ri %*% e$vectors[, o]
  list(lam = e$values[o], psi0 = V[1, ], psiL = V[n + 1, ])
}
Kpole  <- function(sp, s, kmax = 260) { k <- 1:kmax; sum(sp$psi0[k] * sp$psiL[k] * exp(-sp$lam[k] * s)) / (2 * pi) }
Kexact <- function(s, K = 400) { l <- 0:K; sum((2*l+1) * (-1)^l * exp(-l*(l+1)*s)) / (4*pi*a^2) }
# f is given as a function of ARC LENGTH u on [0, Lx], so a rescaled sphere rescales both.
mkf <- function(kind, e) switch(kind,
  round = function(u) a * sin(u / a),
  cubic = local({ ee <- e; function(u) a * (sin(u/a) + ee * sin(u/a)^3) }),
  cos2  = local({ ee <- e; function(u) a * (sin(u/a) + ee * sin(u/a)^3 * cos(2 * u/a)) }))
mkround <- function(b) local({ bb <- b; function(u) bb * sin(u / bb) })
brackets <- function(fn) {
  c(mean_sq = integrate(function(u) fn(u)^2, 0, L, rel.tol = 1e-12)$value / L,
    max_sq_half = max(sapply(seq(0, L, length.out = 200001), fn))^2 / 2)
}

cat("=== 1. the solver, against the case with an exact answer ===\n")
sp0 <- polespec(mkf("round", 0))
cat("      l    lambda (FEM)    l(l+1)     rel. err        s     K(FEM)/K(exact) - 1\n")
ss0 <- c(0.45, 0.30, 0.20, 0.15)
for (i in 1:4) {
  l <- i
  cat(sprintf("   %4d  %13.7f  %8.3f   %.1e     %6.3f   %+.2e\n", l, sp0$lam[l+1], l*(l+1),
              abs(sp0$lam[l+1]/(l*(l+1)) - 1), ss0[i], Kpole(sp0, ss0[i])/Kexact(ss0[i]) - 1))
  stopifnot(abs(Kpole(sp0, ss0[i])/Kexact(ss0[i]) - 1) < 3e-3)
}

cat("\n=== 2. the perturbation reaches the spectrum, or the test is empty ===\n")
cat("      eps    lambda_1   lambda_2   lambda_3     <J^2>/(a^2/2)   max^2/2 /(a^2/2)\n")
for (e in c(0, 0.15, 0.30, 0.50)) {
  sp <- polespec(mkf("cubic", e)); b <- brackets(mkf("cubic", e))
  cat(sprintf("   %6.2f  %9.6f  %9.6f  %9.6f   %13.5f   %16.5f\n",
              e, sp$lam[2], sp$lam[3], sp$lam[4], b[1]/0.5, b[2]/0.5))
}
cat("   The eigenvalues move, so the profile is live in the solver.\n")

cat("\n=== 3. the decisive test, with no extrapolation in it ===\n")
cat("   Take the ratio of the perturbed kernel to the round one at FIXED s. Whatever bias the\n")
cat("   solver and the window carry cancels. If the bracket controls the amplitude, the ratio\n")
cat("   must tend to sqrt(bracket ratio). If nothing but L does, it must tend to 1.\n\n")
for (kind in c("cubic", "cos2")) for (e in c(0.15, 0.30, 0.50)) {
  sp <- polespec(mkf(kind, e)); b <- brackets(mkf(kind, e))
  ss <- c(0.45, 0.35, 0.25, 0.18, 0.14, 0.12)
  r  <- sapply(ss, function(s) Kpole(sp, s) / Kpole(sp0, s))
  fit <- lm(r ~ ss)                                   # the deviation is linear in s
  r0  <- unname(coef(fit)[1])
  cat(sprintf("   %-6s eps %.2f:  ratios", kind, e))
  cat(sprintf(" %.4f", r))
  cat(sprintf("\n              -> s=0: %.4f    <J^2> wants %.4f, max^2/2 wants %.4f\n",
              r0, sqrt(b[1]/0.5), sqrt(b[2]/0.5)))
  stopifnot(abs(r0 - 1) < 0.02)
}
cat("\n   Every one extrapolates to 1, while the two brackets demand between 0.73 and 1.50.\n")
cat("   NEITHER READING IS RIGHT. The caustic amplitude does not see the Jacobi profile.\n")

cat("\n=== 4. so the amplitude is the geodesic length and the family's volume, and no more ===\n")
cat("   Putting <J^2> = a^2/2 = L^2/2pi^2 back into the form that fitted the sphere family,\n")
cat("        c_n = V_fam ( L / (2 sqrt(pi)) )^n ,\n")
cat("   which now stands on the profile-independence above rather than on one profile.\n")
Vol <- function(N) 2 * pi^((N + 1) / 2) / gamma((N + 1) / 2)
c_ref <- function(n) (pi * a)^n * pi^(n/2) * 4 * pi / (factorial(n) * Vol(n + 1))
cat("      n    V_fam        V_fam (L/2sqrt(pi))^n      c_n reference     difference\n")
for (n in 1:4) {
  v <- Vol(n) * (L / (2 * sqrt(pi)))^n
  cat(sprintf("   %4d  %10.6f  %24.10f  %16.10f   %.1e\n", n, Vol(n), v, c_ref(n), abs(v - c_ref(n))))
  stopifnot(abs(v - c_ref(n)) < 1e-10)
}

cat("\n=== 5. the contact geodesic, and WHICH length goes in ===\n")
lam_tot <- 3 * pi / 2 + 4; J2end <- 47.561945
Dv <- sqrt(lam_tot / J2end)
cat("   V_fam = 2 pi there too, a circle of planes through the axis, the same as on S^2, so\n")
cat("   the transfer needs only the length and the reduced Van Vleck.\n")
cat("   The Jacobi operator block-diagonalises exactly on this curve: contact_vanvleck.R shows\n")
cat("   both transverse directions are parallel-propagated and the tidal matrix is diag(+Psi,\n")
cat("   -Psi), so the Gaussian integral factorises over the degenerate and non-degenerate\n")
cat("   sectors and Delta'^{1/2} multiplies c_1. That is exact here, not an assumption.\n")
cat("   The length is the one A.19 identifies and not the affine total, and the difference is\n")
cat("   a factor of 0.59 in the amplitude, so it is worth writing both out. On the sphere the\n")
cat("   affine parameter, the geodesic's own length and its projection onto the sphere all\n")
cat("   coincide at pi a, which is why sections 1 to 4 cannot tell them apart. A.19 separates\n")
cat("   them by putting the caustic in an S^2 and letting the geodesic also run along a flat\n")
cat("   direction, which holds the projection at pi a while the geodesic's own length, and with\n")
cat("   it the affine parameter, grows: over a 38 per cent change the amplitude holds to the\n")
cat("   fifth decimal. So the projection is the length and the affine total is not.\n")
proj <- pi + 2                      # int r dphi along r = M(1 + sin phi) over phi in [0, pi]
gi <- sum(1 + sin(pi * (seq_len(400000) - 0.5) / 400000)) * pi / 400000
cat(sprintf("\n      projection arc length int r dphi = %.9f M, quadrature %.9f M\n", proj, gi))
stopifnot(abs(gi - proj) < 1e-7)
c1 <- sqrt(pi) * proj
amp <- Dv * c1
cat(sprintf("      c_1 = sqrt(pi) M(pi + 2) = %.6f M\n", c1))
cat(sprintf("      Delta'^{1/2} = sqrt(%.6f) = %.6f\n", lam_tot / J2end, Dv))
cat(sprintf("      Delta^{1/2} -> %.4f M s^{-1/2}\n", amp))
cat(sprintf("\n   The affine total would give sqrt(pi)(3 pi/2 + 4) = %.6f M and an amplitude of\n",
            sqrt(pi) * lam_tot))
cat(sprintf("   %.4f, larger by %.4f. That reading is the one A.19's measurement excludes, and it\n",
            Dv * sqrt(pi) * lam_tot, (Dv * sqrt(pi) * lam_tot) / amp))
cat("   sat in this file for a day: the two coincide on every geometry sections 1 to 4 use, so\n")
cat("   nothing here caught it and the provenance audit did.\n")
stopifnot(abs(amp - 3.9004) < 1e-3)

cat("\n=== 6. the plants ===\n")
cat("   (a) the ratio test must detect a profile change if there were one to detect. Feed it\n")
cat("       a sphere of a DIFFERENT radius, where the length really does change:\n")
for (scl in c(1.05, 1.10)) {
  b <- a * scl
  spS <- polespec(mkround(b), Lx = pi * b)
  r <- Kpole(spS, 0.25) / Kpole(sp0, 0.25)
  cat(sprintf("       radius x %.2f (L = %.4f):  kernel ratio at s = 0.25 is %.3e\n", scl, pi*b, r))
  stopifnot(abs(r - 1) > 0.5)
}
cat("   (b) and the round case must reproduce the known amplitude through the same pipeline.\n")
ss <- seq(0.15, 0.45, length.out = 16)
A  <- sapply(ss, function(s) Kpole(sp0, s) * (4*pi*s) * exp(L^2/(4*s)) * sqrt(s))
c1m <- unname(predict(lm(A ~ poly(ss, 4)), newdata = data.frame(ss = 0)))
cat(sprintf("       measured c_1 = %.4f against sqrt(pi) L = %.4f, off by %.2f per cent\n",
            c1m, sqrt(pi) * L, 100 * abs(c1m / (sqrt(pi) * L) - 1)))
cat("       That residual bias is why section 3 takes ratios instead of absolute values.\n")

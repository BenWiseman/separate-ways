#!/usr/bin/env Rscript
# Fork 8's saddle-point route cannot check A.19's replacement rule, because in the limit it works
# in it IS the rule's own input. Worth establishing before anyone spends a session on it.
#
# THE TEMPTATION. The interior's radial modes are known in closed form at zero frequency and in
# the eikonal limit, so the cheap way to the interior sum is to work at large l, saddle-point the
# frequency integral about k = 0, and read the answer off. Fork 8 records that route as untried
# and REORIENT warns it is where a sign gets lost. The deeper problem is not the sign.
#
# WHY IT IS CIRCULAR. The WKB limit of a propagator IS the leading Hadamard form: amplitude the
# square root of the Van Vleck determinant, phase the world function. That is geometric optics and
# it is what the parametrix is. A.18 exists because the parametrix FAILS at a caustic, where the
# Van Vleck determinant diverges, and A.19's rule is the replacement for that divergence. A sum
# taken in the eikonal limit reproduces the object that diverges; it cannot produce the finite
# thing meant to replace it. A.18's geometry was chosen for the opposite property: its reduced
# two-dimensional propagator is EXACT per mode, a Bessel function and not a WKB amplitude, so the
# sum sees the caustic without a Van Vleck factor ever appearing.
#
# This file shows the divergence explicitly along the interior's own contact geodesic, so the
# circularity is a measured statement rather than an argument.

TOL <- 1e-9
fail <- 0
note <- function(ok, what) { if (!ok) { cat("   *** FAILED:", what, "\n"); fail <<- fail + 1 } }

lam_tot <- 3*pi/2 + 4
Dp      <- 0.183180                      # the companion's Delta -> Dp lam_tot/(lam_tot - lambda)
Delta   <- function(lam) Dp*lam_tot/(lam_tot - lam)

cat("=== 1. the eikonal amplitude on the contact geodesic diverges at the caustic ===\n")
cat("      lambda_tot - lambda      Delta        Delta^{1/2}      x (lam_tot-lam)^{1/2}\n")
prod <- c()
for (d in c(1e-1, 1e-2, 1e-3, 1e-4, 1e-6)) {
  lam <- lam_tot - d; D <- Delta(lam)
  cat(sprintf("   %18.0e %14.4f %15.4f %22.6f\n", d, D, sqrt(D), sqrt(D)*sqrt(d)))
  prod <- c(prod, sqrt(D)*sqrt(d))
}
cat(sprintf("   the last column is flat to %.1e, so the eikonal amplitude goes exactly as\n",
            max(abs(diff(prod)))/mean(prod)))
cat("   (lambda_tot - lambda)^{-1/2} and has no finite limit at the caustic.\n")
note(max(abs(diff(prod)))/mean(prod) < 1e-9, "the divergence is a clean inverse square root")
note(abs(prod[1] - sqrt(Dp*lam_tot)) < 1e-9, "with coefficient sqrt(Delta' lambda_tot)")

cat("\n=== 2. and A.19's rule is the finite thing that replaces it ===\n")
cat("   The rule reads Delta^{1/2} -> Delta'^{1/2} sqrt(pi) L s^{-1/2}, with L the projection's\n")
cat("   arc length. Both sides carry an inverse square root, of the DISTANCE TO THE CAUSTIC on\n")
cat("   the left and of the OFFSET s on the right, and the replacement is what relates them.\n")
amp <- sqrt(Dp)*sqrt(pi)*(pi + 2)
cat(sprintf("      Delta'^{1/2} sqrt(pi) M(pi+2) = %.4f M, the target fork 9 must hit\n", amp))
note(abs(amp - 3.9004) < 1e-3, "the rule's value is the calibration target")
cat("   A sum taken in the eikonal limit returns the LEFT side, which is the divergence. It\n")
cat("   cannot return the right side, because the right side is what replaces it.\n")

cat("\n=== 3. the plant: a geometry with no caustic must show no divergence ===\n")
cat("   Off the caustic the Van Vleck factor is finite and the eikonal amplitude is the answer,\n")
cat("   which is why the parametrix is right everywhere else. Take the same form with the pole\n")
cat("   moved off the path, lam_pole = 2 lam_tot:\n")
Dn <- function(lam) Dp*lam_tot/(2*lam_tot - lam)
pn <- sapply(c(1e-1, 1e-3, 1e-6), function(d) sqrt(Dn(lam_tot - d))*sqrt(d))
cat(sprintf("      the same column: %.4f, %.4f, %.4f, which collapses rather than staying flat\n",
            pn[1], pn[2], pn[3]))
note(max(abs(diff(pn))) > 0.01, "plant: with no caustic on the path there is no inverse square root")

cat("\n=== 4. so what fork 9 needs, named ===\n")
cat("   Modes that are accurate AT the caustic, not merely at large l. Two ways in:\n")
cat("   (a) exact radial modes per (k, l), numerically, which is what A.18 had analytically and\n")
cat("       is the expensive route REORIENT records;\n")
cat("   (b) UNIFORM asymptotics rather than eikonal, the standard treatment of a fold caustic,\n")
cat("       where the divergent WKB amplitude is replaced by an integral over the degenerate\n")
cat("       family. A.19's rule already has that shape: L is the family's extent and sqrt(pi)\n")
cat("       s^{-1/2} is the Gaussian in the one non-degenerate direction. So (b) would be a\n")
cat("       second derivation of the rule rather than an independent check of it, and only (a)\n")
cat("       is a check.\n")
cat("   That is worth knowing before starting: the cheap route reproduces the rule by\n")
cat("   construction, and the check has to be the expensive one.\n")

if (fail == 0) cat("\n   all checks passed\n") else {
  cat(sprintf("\n   %d CHECK(S) FAILED\n", fail)); quit(status = 1)
}

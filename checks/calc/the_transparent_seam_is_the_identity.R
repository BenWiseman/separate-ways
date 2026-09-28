#!/usr/bin/env Rscript
# The transparent seam is not one choice out of a family. It is the family's identity element,
# which is to say the absence of a seam.
#
# WHY THIS MATTERS. The seam coefficient carries the paper's sharpest near-term test, an exactly
# Kerr ringdown decidable at about twice present exposure, and it is the one place the algebra has
# not paid: three routes to deriving kappa = 1 are closed and what selects it is a soft
# variational maximum. That reads as "a value was picked out of a one-parameter family", which is
# the weakest possible framing and is also not what happened.
#
# A.11 already establishes the group structure and uses it only for a route that fails. Writing
# kappa = e^u the family's transfer matrix is the boost (cosh u, sinh u; sinh u, cosh u), so seams
# in series add rapidities and kappa -> 1/kappa is the group inverse. The consequence it does not
# draw is the obvious one: u = 0, which is kappa = 1, is the IDENTITY of that group. Adopting the
# transparent seam is adopting no boundary term at all, which is the same minimal choice the
# implementation makes everywhere else, and it puts the burden where it belongs. A reflecting seam
# is an addition, and nothing in the construction calls for one.

TOL <- 1e-12
fail <- 0
note <- function(ok, what) { if (!ok) { cat("   *** FAILED:", what, "\n"); fail <<- fail + 1 } }

Tm <- function(u) matrix(c(cosh(u), sinh(u), sinh(u), cosh(u)), 2, 2)
rk <- function(k) (1 - k^2)/(1 + k^2)
tk <- function(k) 2*k/(1 + k^2)

cat("=== 1. the family is a one-parameter group, and composition adds rapidities ===\n")
cat("      u1        u2       max |T(u1) T(u2) - T(u1+u2)|\n")
for (p in list(c(0.3, 0.7), c(-1.2, 0.4), c(2.0, 2.0), c(0.0, 1.1))) {
  d <- max(abs(Tm(p[1]) %*% Tm(p[2]) - Tm(p[1] + p[2])))
  cat(sprintf("   %8.2f %9.2f %28.2e\n", p[1], p[2], d))
  note(d < TOL, "seams in series add rapidities")
}

cat("\n=== 2. and its identity is u = 0, which is kappa = 1 ===\n")
cat(sprintf("      T(0) - I:              %.2e\n", max(abs(Tm(0) - diag(2)))))
cat(sprintf("      kappa = e^0:           %.6f\n", exp(0)))
cat(sprintf("      r at kappa = 1:        %.2e      t: %.6f\n", rk(1), tk(1)))
note(max(abs(Tm(0) - diag(2))) < TOL, "T(0) is the identity matrix")
note(abs(rk(1)) < TOL && abs(tk(1) - 1) < TOL, "kappa = 1 reflects nothing and transmits everything")
cat("   So the transparent point is the element that composes with anything and changes it not at\n")
cat("   all, which is what 'no seam' means as an operation. Every other kappa is an addition.\n")

cat("\n=== 3. the plant: no other member is an identity, and the test can see that ===\n")
cat("      kappa      u       max |T(u) T(0.5) - T(0.5)|      is it an identity?\n")
for (k in c(0.25, 0.5, 1, 2, 4)) {
  u <- log(k); d <- max(abs(Tm(u) %*% Tm(0.5) - Tm(0.5)))
  cat(sprintf("   %8.3f %8.4f %26.4f %22s\n", k, u, d,
              ifelse(d < TOL, "yes", "no")))
  note((abs(k - 1) < TOL) == (d < TOL), sprintf("only kappa = 1 acts as the identity (k = %g)", k))
}

cat("\n=== 4. what that changes about the claim ===\n")
cat("   Not the value, which is the same, and not the three closed routes, which stay closed. What\n")
cat("   changes is where the burden sits. The implementation's other choices are all of the form\n")
cat("   'add no structure': no horizon reflectivity, no mode mixing, no spectral action. The seam\n")
cat("   is the same choice and was being described as though it were different, a value selected\n")
cat("   from a continuum. It is the continuum's identity. A reflecting seam is an added boundary\n")
cat("   term with its own matching law, and the question a reader should ask is not why this value\n")
cat("   but what would require any.\n")
cat("   What this does NOT do: it does not derive kappa = 1. A theory may add a boundary term, and\n")
cat("   identity elements are not privileged by anything except minimality. The soft variational\n")
cat("   maximum and the corner term's uniqueness among scale-free seams are still what support it.\n")

if (fail == 0) cat("\n   all checks passed\n") else {
  cat(sprintf("\n   %d CHECK(S) FAILED\n", fail)); quit(status = 1)
}

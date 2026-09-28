#!/usr/bin/env Rscript
# world_function_zero.R -- A.15 conjectures that the fold censors its own closed causal
# curves, on the ground that the angular deficit has a simple zero at r = M. That is the
# right kind of evidence but not the quantity an image sum needs. What the image term of
# a free field's two-point function sees is the WORLD FUNCTION sigma(x, Theta x), and
# the question is whether THAT has a simple zero. This computes it.
#
# The connecting curve lies at X = 0, so Schwarzschild t is constant along it and the
# conserved energy vanishes: E = (1 - 2M/r) dt/dtau = 0. The E = 0 timelike geodesics
# inside the horizon therefore obey
#     (dr/dtau)^2 = (2M/r - 1)(1 + L^2/r^2),      dphi/dtau = L/r^2,
# so over the two interior legs, from r up to the horizon and back,
#     Dphi(r,L) = 2 int_r^{2M} (L/r'^2) dr' / sqrt((2M/r'-1)(1 + L^2/r'^2)),
#     tau(r,L)  = 2 int_r^{2M}          dr' / sqrt((2M/r'-1)(1 + L^2/r'^2)).
# Fix L by Dphi = pi, which is what the antipodal map demands, and read off tau.
# As L -> infinity the first integrand tends to 1/sqrt(r'(2M-r')) and Dphi tends to the
# null bound 2 pi - 4 arcsin sqrt(r/2M), so the L -> infinity limit IS the r = M case.
#
# sigma = -tau^2/2 for a timelike pair, so a square-root vanishing of tau is a SIMPLE
# zero of sigma, which is what a power-law divergence of the image sum requires.

M <- 1; rh <- 2*M
# substitute r' = rh - u^2 to kill the inverse-square-root endpoint exactly
integ <- function(r, L, kind) {
  f <- function(u) { x <- rh - u*u
    base <- (rh/x - 1)*(1 + L^2/x^2)
    num <- if (kind == "phi") L/x^2 else 1
    2*u*num/sqrt(pmax(1e-300, base)*1) * (1/u) }   # dr' = 2u du, and 1/u from sqrt(2M/x-1)
  # (2M/x - 1) = (rh-x)/x = u^2/x, so sqrt(base) = (u/sqrt(x)) sqrt(1+L^2/x^2)
  g <- function(u) { x <- rh - u*u
    num <- if (kind == "phi") L/x^2 else 1
    2*num*sqrt(x)/sqrt(1 + L^2/x^2) }
  2*integrate(g, 0, sqrt(rh - r), rel.tol=1e-12, subdivisions=3000L)$value
}
dphi <- function(r, L) integ(r, L, "phi")
tau  <- function(r, L) integ(r, L, "tau")

cat("=== 1. the null limit reproduces the known bound ===\n")
cat("      r/M      Dphi at L=1e7      2pi - 4 asin sqrt(r/2M)      diff\n")
for (x in c(0.2, 0.5, 0.8, 1.0)) {
  r <- x*M; a <- dphi(r, 1e7); b <- 2*pi - 4*asin(sqrt(r/rh))
  cat(sprintf("   %6.2f       %10.7f              %10.7f       %.1e\n", x, a, b, a-b))
  stopifnot(abs(a-b) < 1e-5)
}
cat("   So the large-L family limits onto the null curve, as it must.\n")

cat("\n=== 2. the angular momentum that delivers exactly pi, and the proper time ===\n")
L_for_pi <- function(r) uniroot(function(l) dphi(r, exp(l)) - pi,
                                c(-4, 20), tol=1e-12)$root
cat("      r/M        L          tau (M)       tau^2/(M-r)     |sigma| = tau^2/2\n")
prev <- NULL
for (x in c(0.10, 0.30, 0.50, 0.70, 0.90, 0.97, 0.99, 0.997, 0.999)) {
  r <- x*M; L <- exp(L_for_pi(r)); tt <- tau(r, L)
  cat(sprintf("   %6.3f  %9.3f    %10.6f    %12.6f    %12.6f\n",
              x, L, tt, tt^2/(M-r), tt^2/2))
  prev <- c(x, tt)
}
cat("\n   tau^2/(M-r) tends to a finite nonzero constant, so tau ~ sqrt(M-r) and\n")
cat("   sigma = -tau^2/2 has a SIMPLE zero at r = M.\n")
k <- sapply(c(0.99, 0.997, 0.999, 0.9997), function(x) {
  r <- x*M; L <- exp(L_for_pi(r)); tau(r,L)^2/(M-r) })
cat(sprintf("   limiting slope from the last four rows: %s\n", paste(sprintf("%.4f", k), collapse=" ")))
stopifnot(all(is.finite(k)), max(k)/min(k) < 1.25, min(k) > 0.1)

cat("\n=== 3. the deepest loop, for the record ===\n")
cat("      r0/M        proper length (M)\n")
for (r0 in c(1e-2, 1e-3, 1e-4, 1e-5, 1e-6)*M) {
  L0 <- exp(L_for_pi(r0))
  cat(sprintf("   %9.0e        %12.6f\n", r0/M, tau(r0, L0)))
}
r0 <- 1e-6*M; L0 <- exp(L_for_pi(r0))
cat(sprintf("   the constant-winding ansatz sqrt(3) pi M = %.6f M is only a lower bound\n", sqrt(3)*pi))

cat("\n=== 4. PLANTED FAILURES ===\n")
b1 <- sapply(c(0.99, 0.999), function(x) { r <- x*M; L <- exp(L_for_pi(r)); tau(r,L)^2/(M-r)^2 })
cat(sprintf("   (a) testing for a DOUBLE zero, tau^2/(M-r)^2, gives %.1f then %.1f: diverging,\n", b1[1], b1[2]))
cat("       so the zero is simple and not double. The test can tell them apart.\n")
stopifnot(b1[2]/b1[1] > 5)
b2 <- sapply(c(0.99, 0.999), function(x) { r <- x*M; L <- exp(L_for_pi(r)); tau(r,L)/(M-r) })
cat(sprintf("   (b) testing tau itself against (M-r) gives %.1f then %.1f: also diverging,\n", b2[1], b2[2]))
cat("       so tau does not vanish linearly either. Only tau^2 ~ (M-r) survives.\n")
stopifnot(b2[2]/b2[1] > 2)

cat(sprintf("
=== flatly ===

  The world function has a simple zero. Fixing the angular momentum so the two interior
  legs deliver exactly pi, the proper time of the connecting geodesic vanishes as
  sqrt(M - r), so sigma = -tau^2/2 vanishes LINEARLY at r = M. A double zero and a
  linear vanishing of tau are both excluded by the same numbers.

  That is what the conjecture needed. The image term of a free field's two-point
  function on the quotient goes as 1/sigma, so it diverges as (M-r)^{-1} on approach to
  the contact boundary, and the stress tensor, built from two derivatives of it, as a
  higher inverse power. We do not fix the power here: that needs the Hadamard
  coefficients of the actual mode functions and not just the world function, and the
  sign, which decides whether the divergence closes the region or merely marks it, needs
  the same. What is now computed is the input, and it is the input that was in doubt.\n"))

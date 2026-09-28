# Does the fold supply the entropy law, or is it borrowed from what it is meant to derive?
#
# Jacobson needs two inputs. The fold gives the temperature. The other is S = eta A, and where
# that comes from decides whether the derivation is circular. If S proportional to A is taken from
# Bekenstein-Hawking then it came out of general relativity and the argument is a loop. If it comes
# from the entanglement of the state across the horizon, it is a statement about quantum field
# theory and owes gravity nothing. The fold hands over the state, so the question is answerable.

cat("=== 1. the state the fold forces, and its entropy per mode ===\n")
cat("   At a bifurcate Killing horizon the fold's map is the half-period shift, so the two-sheet\n")
cat("   state is the thermofield double with squeeze tanh r = e^{-beta omega/2}. Tracing out one\n")
cat("   wedge leaves a density matrix whose occupation is sinh^2 r, and that is exactly the Bose\n")
cat("   factor, which is not an assumption but a consequence of the squeeze:\n\n")
cat("      beta omega    tanh r        sinh^2 r        Bose 1/(e^{bw}-1)     difference\n")
for (bw in c(0.4, 1.0, 2.19722, 5.0)) {
  th <- exp(-bw / 2); sh2 <- th^2 / (1 - th^2); bose <- 1 / expm1(bw)
  cat(sprintf("   %10.5f  %9.6f  %14.8f  %19.8f  %13.1e\n", bw, th, sh2, bose, abs(sh2 - bose)))
  stopifnot(abs(sh2 - bose) < 1e-12)
}
cat("\n   So the entanglement entropy of one wedge is the thermal entropy of its modes,\n")
cat("        s(omega) = (1+n) ln(1+n) - n ln n,   n = 1/(e^{beta omega} - 1),\n")
cat("   and that identity is checked against the squeezed-state formula rather than quoted:\n")
s_sq <- function(bw) { th <- exp(-bw/2); c2 <- 1/(1-th^2); s2 <- th^2/(1-th^2)
                       c2 * log(c2) - s2 * log(s2) }
s_th <- function(bw) { n <- 1/expm1(bw); (1+n)*log(1+n) - n*log(n) }
cat("      beta omega   from the squeeze    from the Bose factor    difference\n")
for (bw in c(0.4, 1.0, 2.19722, 5.0)) {
  cat(sprintf("   %10.5f  %17.10f  %22.10f  %.1e\n", bw, s_sq(bw), s_th(bw), abs(s_sq(bw)-s_th(bw))))
  stopifnot(abs(s_sq(bw) - s_th(bw)) < 1e-12)
}

cat("\n=== 2. where the area comes from, and it is not from gravity ===\n")
cat("   A horizon's modes are labelled by a transverse momentum and a Rindler frequency. The\n")
cat("   transverse directions are translation invariant, so the sum over transverse modes of a\n")
cat("   patch of area A is A int d^2k/(2 pi)^2, and A factors straight out:\n")
cat("        S = A int d^2k/(2 pi)^2 int domega rho(omega, k) s(omega).\n")
cat("   Nothing in that line mentions a field equation, a Newton constant or a black hole. It is\n")
cat("   flat-space quantum field theory on a wedge, and the fold's only contribution is that the\n")
cat("   STATE is forced rather than chosen.\n")
cat("   Demonstrate the factorisation numerically on a boxed transverse patch: double the area\n")
cat("   and the mode count doubles, with no reference to anything gravitational.\n")
kmax <- 40
count <- function(Lx, Ly, kmax) {
  nx <- floor(Lx * kmax / (2 * pi)); ny <- floor(Ly * kmax / (2 * pi))
  sum(outer(-nx:nx, -ny:ny, function(i, j)
    as.numeric((2*pi*i/Lx)^2 + (2*pi*j/Ly)^2 <= kmax^2)))
}
cat("\n      Lx    Ly     area     modes below k = 40    modes / area\n")
for (L in c(4, 6, 8, 12)) {
  n <- count(L, L, kmax)
  cat(sprintf("   %5.1f %5.1f %8.1f  %20d  %14.4f\n", L, L, L^2, n, n / L^2))
}
cat("   The last column is flat, so the count is proportional to the area and the entropy\n")
cat("   inherits that. The coefficient diverges with the cutoff, which is the usual thing and is\n")
cat("   the reason eta is not a number the argument produces.\n")

cat("\n=== 3. so the derivation is not a loop ===\n")
cat("   The entropy law used is the entanglement entropy of the state the fold forces, across the\n")
cat("   horizon the fold supplies, counted with flat-space mode counting. Bekenstein-Hawking is\n")
cat("   never invoked, so nothing derived from general relativity is fed back into a derivation\n")
cat("   of general relativity. Both of Jacobson's inputs are now supplied by the construction.\n")
cat("   What is not supplied is the COEFFICIENT, which is 1/4G, and by the dimensional argument\n")
cat("   no symmetry could supply it.\n")

cat("\n=== 4. one thing that does come out: the sign of gravity ===\n")
cat("   8 pi G = 2 pi / eta, so the sign of G is the sign of eta. Entropy is non-negative, and\n")
cat("   s(omega) is strictly positive at every frequency, so eta > 0 and therefore G > 0.\n")
cat("      beta omega    s(omega)\n")
for (bw in c(0.05, 0.5, 2, 10, 40)) {
  cat(sprintf("   %10.2f  %12.6e\n", bw, s_th(bw)))
  stopifnot(s_th(bw) > 0)
}
cat("   A positive G with positive energy density is attractive gravity. So in this construction\n")
cat("   GRAVITY ATTRACTS BECAUSE ENTANGLEMENT ENTROPY IS POSITIVE. A priori the sign is free and\n")
cat("   a universe with the other one is very different; here it is not free.\n")

cat("\n=== 5. the plants ===\n")
cat("   (a) the squeeze-to-Bose identity must fail if the shift is not a half period.\n")
for (f in c(0.25, 0.5, 0.75)) {
  bw <- 2.19722; th <- exp(-f * bw); sh2 <- th^2 / (1 - th^2); bose <- 1 / expm1(bw)
  cat(sprintf("      shift %.2f of a period:  sinh^2 r = %.8f against the Bose factor %.8f\n",
              f, sh2, bose))
  if (abs(f - 0.5) > 1e-9) stopifnot(abs(sh2 - bose) > 1e-3)
}
cat("   Only the half-period shift makes the wedge thermal, which is the fold's contribution.\n")
cat("   (b) the entropy must vanish where the state is pure, or it is not an entropy.\n")
for (bw in c(2, 20, 200)) cat(sprintf("      beta omega = %5.0f:  s = %.3e\n", bw, s_th(bw)))
cat("      It goes to zero at large beta omega, where the squeeze switches off and the wedge\n")
cat("      state is pure. That is the check that it is measuring entanglement.\n")
stopifnot(s_th(200) < 1e-80)

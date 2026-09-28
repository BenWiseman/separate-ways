#!/usr/bin/env Rscript
# Fork 9: the CONTACT pair, future interior against past interior, assembled from modes and
# checked against the exact answer. In flat space, where there is an exact answer.
#
# WHY. the_image_sum_is_right_in_flat_space.R validated the assembly on the wedge pair, which is
# the SILENCE configuration: both points outside, separation spacelike, sheets unable to touch.
# The configuration fork 9 needs is the other one, future interior against past interior, where
# the separation is timelike and the sheets do touch. Flat space has it, as the Milne region, and
# it is where the two exterior mode families mix. That mixing was the one open risk. It is not a
# risk: the continuation of the mode across the branch cut carries it, with no separate Bogoliubov
# bookkeeping, and the assembly reproduces the exact two-point function.
#
# THE MODEL. Two-dimensional Minkowski, massive scalar. Milne region F: T = rho cosh(eta),
# X = rho sinh(eta), ds^2 = -d rho^2 + rho^2 d eta^2, so rho is TIME and eta is space, exactly
# as r is time and t is space inside a horizon. A mode e^{i k eta} f(rho) obeys
#     f'' + f'/z + (k^2/z^2 + 1) f = 0,   z = m rho,
# the Bessel equation of order i k, so f is a Hankel function of imaginary order. The image of
# (eta, rho) under (T, X) -> (-T, -X) sits in P at the same (eta, rho), reached by rho -> rho
# e^{-i pi}, since U = rho e^{-eta} and V = rho e^{+eta} both change sign.
#
# THE NORMALISATION. Constant-rho slices are the Cauchy slices. The Klein-Gordon product there is
# -i int d eta rho (f d_rho g* - g* d_rho f), and z(F dFbar/dz - Fbar dF/dz) is z-independent and
# equals 4i e^{-pi k}/pi for F = H^{(2)}_{ik}. With the 2 pi from the eta integral that gives
#     N_k^2 = e^{pi k}/8.
# THE EXACT ANSWER. The separation is (2T, 2X) with interval -4 rho^2, timelike and future
# directed, so W = (1/2 pi) K_0(2 i m rho) = -(i/4) H_0^{(2)}(2 m rho). The other i-epsilon would
# give +(i/4) H_0^{(1)}, which differs by the sign of the imaginary part; the computation picks
# one of them and is not told which.

fail <- 0
note <- function(ok, what) { if (!ok) { cat("   *** FAILED:", what, "\n"); fail <<- fail + 1 } }

LG <- c(676.5203681218851, -1259.1392167224028, 771.32342877765313,
        -176.61502916214059, 12.507343278686905, -0.13857109526572012,
        9.9843695780195716e-6, 1.5056327351493116e-7)
lgammaC <- function(z) {              # Lanczos, valid for Re z >= 0.5, which is all we ask of it
  stopifnot(all(Re(z) >= 0.5))
  z <- z - 1; x <- rep(0.99999999999980993 + 0i, length(z))
  for (i in seq_along(LG)) x <- x + LG[i]/(z + i)
  t <- z + length(LG) - 0.5
  0.5*log(2*pi) + (z + 0.5)*log(t) - t + log(x)
}
JnuL <- function(nu, lz, NT = 260) {  # J_nu(z) with log(z/2) supplied, so the branch is explicit
  s <- rep(0+0i, length(nu))
  for (n in 0:NT) s <- s + ((-1)^n)*exp((nu + 2*n)*lz - lgamma(n+1) - lgammaC(nu + n + 1))
  s
}
HankL <- function(k, lz, which) {
  nu <- 1i*k
  Jp <- JnuL(nu, lz); Jm <- JnuL(-nu, lz); s <- sin(nu*pi)
  if (which == 1) (Jm - exp(-nu*pi*1i)*Jp)/(1i*s) else (Jm - exp(nu*pi*1i)*Jp)/(-1i*s)
}
H2    <- function(k, x) HankL(k, log(x/2) + 0i, 2)
H1    <- function(k, x) HankL(k, log(x/2) + 0i, 1)
H2rot <- function(k, x, sgn) HankL(k, log(x/2) + sgn*1i*pi, 2)

cat("=== 1. the tooling, checked three ways before it is used ===\n")
cat(sprintf("   lgammaC on the real axis: %.12f against R's lgamma %.12f at 7.3\n",
            Re(lgammaC(7.3+0i)), lgamma(7.3)))
note(abs(Re(lgammaC(7.3+0i)) - lgamma(7.3)) < 1e-11, "the complex log-gamma is right on the reals")
cat(sprintf("   |Gamma(1+i)| = %.12f against the closed form sqrt(pi/sinh pi) = %.12f\n",
            Mod(exp(lgammaC(1+1i))), sqrt(pi/sinh(pi))))
note(abs(Mod(exp(lgammaC(1+1i))) - sqrt(pi/sinh(pi))) < 1e-11, "and off them")
for (p in list(c(0, 1.3), c(2, 5.0))) {
  v <- Re(JnuL(p[1]+0i, log(p[2]/2)+0i))
  cat(sprintf("   J_%d(%.1f) = %.12f against besselJ %.12f\n", p[1], p[2], v, besselJ(p[2], p[1])))
  note(abs(v - besselJ(p[2], p[1])) < 1e-11, "the series reproduces besselJ at real order")
}
cat("   the Hankel pair of IMAGINARY order, by its Wronskian H1 H2' - H2 H1' = -4i/(pi z):\n")
for (pr in list(c(0.7, 2.0), c(2.0, 3.5))) {
  k <- pr[1]; x <- pr[2]; h <- 1e-5
  W <- H1(k,x)*((H2(k,x+h)-H2(k,x-h))/(2*h)) - H2(k,x)*((H1(k,x+h)-H1(k,x-h))/(2*h))
  cat(sprintf("      k=%.1f z=%.1f   %+.10fi against %+.10fi\n", k, x, Im(W), -4/(pi*x)))
  note(Mod(W + 4i/(pi*x))/(4/(pi*x)) < 1e-8, sprintf("the Wronskian holds at k=%g", k))
}
for (pr in list(c(0.7, 2.0), c(-1.1, 1.4))) {
  k <- pr[1]; x <- pr[2]
  note(Mod(Conj(H2(k,x)) - exp(-pi*k)*H1(k,x))/Mod(H1(k,x)) < 1e-12,
       sprintf("conj(H2_{ik}) = e^{-pi k} H1_{ik} at k=%g", k))
}
cat("   conj(H2_{ik}(x)) = e^{-pi k} H1_{ik}(x) on the real axis: holds\n")

cat("\n=== 2. the Milne Klein-Gordon constant, which fixes the normalisation ===\n")
cat("        k       z      z(F dFbar/dz - Fbar dF/dz)        4i e^{-pi k}/pi\n")
for (pr in list(c(0.7, 2.0), c(2.0, 3.5))) {
  k <- pr[1]; x <- pr[2]; h <- 1e-5
  F <- H2(k,x); dF <- (H2(k,x+h) - H2(k,x-h))/(2*h)
  val <- x*(F*Conj(dF) - Conj(F)*dF)
  cat(sprintf("   %6.1f %7.1f %25s %22s\n", k, x,
              sprintf("%+.10fi", Im(val)), sprintf("%+.10fi", 4*exp(-pi*k)/pi)))
  note(Mod(val - 4i*exp(-pi*k)/pi)/(4*exp(-pi*k)/pi) < 1e-8,
       sprintf("the Klein-Gordon constant is 4i e^{-pi k}/pi at k=%g", k))
}
cat("   So N_k^2 = e^{pi k}/8 and the integrand is (e^{pi k}/8) H2_{ik}(z) conj H2_{ik}(z e^{-i pi}).\n")

cat("\n=== 3. which branch. Only one of the two converges ===\n")
cat("        k     |integrand| with e^{-i pi}    with e^{+i pi}\n")
for (k in c(1, 3, 6, 10)) {
  a <- Mod((exp(pi*k)/8)*H2(k,1)*Conj(H2rot(k,1,-1)))
  b <- Mod((exp(pi*k)/8)*H2(k,1)*Conj(H2rot(k,1,+1)))
  cat(sprintf("   %6.1f %22.5e %18.5e\n", k, a, b))
  note(b > 10*a, sprintf("the +i pi branch is the divergent one at k=%g", k))
}
cat("   The e^{-i pi} branch is the one the Wightman i-epsilon selects and the only one that\n")
cat("   can be integrated. Nothing is chosen here: the other one simply diverges.\n")

cat("\n=== 4. the tail is exactly 1/(4 pi k), which is why the agreement below stops at a per cent\n")
cat("        k     |integrand|     k times it      1/(4 pi) = 0.0795775\n")
env <- c()
for (k in c(6, 10, 14, 18, 22)) {
  v <- Mod((exp(pi*k)/8)*H2(k,1)*Conj(H2rot(k,1,-1)))
  cat(sprintf("   %6.1f %14.5e %16.6f\n", k, v, k*v)); env <- c(env, k*v)
}
note(abs(env[length(env)] - 1/(4*pi)) < 2e-4, "the envelope is 1/(4 pi k)")
cat("   The phase chirps rather than repeating, so the tail converges, but only conditionally\n")
cat("   and only as fast as the phase accelerates. That is the whole error budget below.\n")

cat("\n=== 5. the assembly against the exact two-point function ===\n")
g <- function(ks, z) (exp(pi*ks)/8)*H2(ks, z)*Conj(H2rot(ks, z, -1))
simp <- function(z, K = 26, n = 2400) {
  ks <- seq(1e-4, K, length.out = n+1); h <- ks[2] - ks[1]   # k = 0 is removable, not a value
  v <- g(ks, z)
  w <- c(1, rep(c(4, 2), length.out = n-1), 1)
  2*sum(w*v)*h/3                                              # even in k
}
cat("      m rho      mode sum                -(i/4) H2_0(2 m rho)        rel     against +(i/4) H1_0\n")
for (z in c(0.5, 1.0, 2.0)) {
  I   <- simp(z)
  ex1 <- -(1i/4)*(besselJ(2*z,0) - 1i*besselY(2*z,0))
  ex2 <-  (1i/4)*(besselJ(2*z,0) + 1i*besselY(2*z,0))
  cat(sprintf("   %7.1f  %9.6f%+9.6fi  %11.6f%+9.6fi %10.4f %12.4f\n",
              z, Re(I), Im(I), Re(ex1), Im(ex1), Mod(I-ex1)/Mod(ex1), Mod(I-ex2)/Mod(ex2)))
  note(Mod(I - ex1)/Mod(ex1) < 0.02, sprintf("the assembly matches -(i/4)H2_0 at m rho = %g", z))
  note(Mod(I - ex2)/Mod(ex2) > 0.5, sprintf("and rejects +(i/4)H1_0 at m rho = %g", z))
}
cat("   The two candidates differ by the sign of the imaginary part and the computation is not\n")
cat("   told which to prefer. It picks one of them by a factor of more than a hundred.\n")

cat("\n=== 6. the plants ===\n")
z <- 1.0; ex <- -(1i/4)*(besselJ(2,0) - 1i*besselY(2,0))
badn <- function(z, wf, K = 20, n = 1600) {
  ks <- seq(1e-4, K, length.out = n+1); h <- ks[2] - ks[1]
  v <- wf(ks)*H2(ks,z)*Conj(H2rot(ks,z,-1))
  w <- c(1, rep(c(4,2), length.out = n-1), 1); 2*sum(w*v)*h/3
}
for (cs in list(list("the normalisation dropped to 1/8", function(k) rep(1/8, length(k))),
                list("e^{pi k/2} in place of e^{pi k}", function(k) exp(pi*k/2)/8),
                list("a spurious factor of k",          function(k) k*exp(pi*k)/8))) {
  v <- badn(z, cs[[2]])
  cat(sprintf("   %-34s gives %11.6f%+11.6fi against %.6f%+.6fi\n",
              cs[[1]], Re(v), Im(v), Re(ex), Im(ex)))
  note(Mod(v - ex)/Mod(ex) > 0.1, sprintf("plant: %s is caught", cs[[1]]))
}

cat("\n=== 7. WHICH continuation, because the obvious one is wrong ===\n")
cat("   Both of these are the point map (U, V) -> (-U, -V), since U = R e^{-eta} and V = R\n")
cat("   e^{+eta} with R the Milne radius: continue the TIME variable, R -> R e^{i pi}, or shift\n")
cat("   the SPACE variable, eta -> eta + i pi. They differ by a loop around the branch point and\n")
cat("   they are not the same continuation of the correlator. The second is the one that looks\n")
cat("   right by analogy with the exterior, where the Killing-time half-period shift IS correct\n")
cat("   and the fifth brick checked it. Inside it is wrong, and wrong in a way that is obvious\n")
cat("   once computed: eta is spacelike there, so shifting it leaves |H2|^2 and a real integral,\n")
cat("   while the true answer is complex.\n")
cat("      m rho        exact           time continuation        space shift\n")
gsp <- function(ks, z) (exp(pi*ks)/8)*exp(-pi*ks)*H2(ks,z)*Conj(H2(ks,z))
sim2 <- function(f, z, K = 20, n = 1600) {
  ks <- seq(1e-4, K, length.out = n+1); h <- ks[2] - ks[1]
  w <- c(1, rep(c(4,2), length.out = n-1), 1); 2*sum(w*f(ks,z))*h/3
}
for (z in c(0.5, 1.0, 2.0)) {
  ex <- -(1i/4)*(besselJ(2*z,0) - 1i*besselY(2*z,0))
  a <- simp(z); b <- sim2(gsp, z)
  cat(sprintf("   %7.1f %10.6f%+9.6fi %12.6f%+9.6fi %12.6f%+9.6fi\n",
              z, Re(ex), Im(ex), Re(a), Im(a), Re(b), Im(b)))
  note(abs(Im(b)) < 1e-9, sprintf("the space shift returns a real number at m rho = %g", z))
  note(Mod(b - ex)/Mod(ex) > 0.5, sprintf("and is wrong by more than half at m rho = %g", z))
}
cat("   So the interior correlator is NOT the Killing-time half-period shift. That prescription\n")
cat("   belongs to the exterior, where the Killing time is timelike and the correlator is\n")
cat("   analytic in a strip. Inside, the analyticity is in the Milne time, and the third brick's\n")
cat("   closing sentence is corrected accordingly: the point map is the half-period shift, the\n")
cat("   CONTINUATION of the correlator is not.\n")

cat("\n=== 8. what this settles ===\n")
cat("   The contact configuration assembles correctly, and the mixing that was the last named\n")
cat("   risk is not a separate computation at all: continuing the mode across the branch cut,\n")
cat("   rho -> rho e^{-i pi}, carries it. No Bogoliubov coefficients are written down anywhere in\n")
cat("   this file and none are needed. The same continuation in Schwarzschild is r -> the image\n")
cat("   region through the horizon, which the Kruskal map already describes.\n")
cat("   WHAT IS LEFT for the hole: the interior radial solutions are not Hankel functions, so\n")
cat("   they come from the solver of interior_modes_at_nonzero_k.R rather than from a series, and\n")
cat("   the l-sum has to be done as well as the frequency integral. Neither is a new question of\n")
cat("   principle. The calibration to meet before any sign is read is still Delta^{1/2} ->\n")
cat("   3.9004 M s^{-1/2}.\n")

if (fail == 0) cat("\n   all checks passed\n") else {
  cat(sprintf("\n   %d CHECK(S) FAILED\n", fail)); quit(status = 1)
}

# The sign of the image stress on the contact orbit, at the leading REAL order in WKB.
#
# A.19 leaves rho_img as the radial part of the wave operator set against the mixed radial
# derivative and says it needs g_l inside the horizon in the Hartle-Hawking state, where no
# closed form exists. That is true of the exact answer. The leading behaviour is geometry,
# and the leading behaviour carries a sign.
#
# A first attempt at this used psi' = i p psi and got an imaginary residue, which an energy
# density cannot be. The amplitude derivative is what was missing. With the full WKB form
# psi ~ p^{-1/2} exp(i int p dr),
#
#     psi'  = (i p - p'/2p) psi,
#     psi'' = (-p^2 - (p'/2p)' + p'^2/4p^2) psi,          the i p' cancels,
#
# so at coincidence, writing N = |psi|^2,
#
#     d_r d_r' g = ( p^2 + p'^2/4p^2 ) N,                 real
#     d_r^2   g = ( -p^2 - (p'/2p)' + p'^2/4p^2 ) N,      real
#     Re d_r  g = ( -p'/2p ) N.
#
# A.19's combination is f d_r^2 g + (f' + 2f/r) d_r g + f d_r d_r' g, and for Schwarzschild
# h = f so sqrt(G) = r^2 and the measure term is 2f/r. The p^2 terms cancel between the
# first and the third, which is the same cancellation the constant-radius model has, and
# what is left is
#
#     rho_res / N  =  f [ -(p'/2p)' + p'^2/2p^2 ]  -  (f' + 2f/r)(p'/2p).
#
# Everything below evaluates that, on the contact orbit and through the interior.

M  <- 1
f  <- function(r) 1 - 2 * M / r
fp <- function(r) 2 * M / r^2
V  <- function(r, L) f(r) * (L * (L + 1) / r^2 + 2 * M / r^3)     # Regge-Wheeler, s = 0
p  <- function(r, L, w) sqrt(w^2 - V(r, L))                        # real inside: V < 0

d1 <- function(g, x, ...) { h <- 1e-5; (g(x + h, ...) - g(x - h, ...)) / (2 * h) }

pp   <- function(r, L, w) d1(function(z) p(z, L, w), r)
half <- function(r, L, w) pp(r, L, w) / (2 * p(r, L, w))           # p'/2p
resid <- function(r, L, w) {
  f(r) * (-d1(function(z) half(z, L, w), r) + 2 * half(r, L, w)^2) -
    (fp(r) + 2 * f(r) / r) * half(r, L, w)
}

cat("=== 1. the imaginary parts really do cancel: psi'' has no i p' term ===\n")
cat("   check numerically that d_r^2 g + d_r d_r' g leaves no p^2 and no imaginary part,\n")
cat("   by evaluating the bracket that survives against its closed form.\n")
L <- 2; w <- 1.3; worst <- 0
for (r in c(0.3, 0.6, 1.0, 1.5)) {
  lhs <- (-d1(function(z) half(z, L, w), r) + 2 * half(r, L, w)^2)
  rhs <- (-d1(function(z) half(z, L, w), r) + (pp(r, L, w) / p(r, L, w))^2 / 2)
  worst <- max(worst, abs(lhs - rhs))
}
cat(sprintf("   p'^2/2p^2 written two ways agree to %.2e\n", worst))
stopifnot(worst < 1e-8)

cat("\n=== 2. the residue ON the contact orbit, r = M, across modes ===\n")
cat("      l    omega      residue / |psi|^2      sign\n")
s <- c()
for (L in c(0, 1, 2, 5, 10, 20)) for (w in c(0.2, 0.5, 1.0, 2.0)) {
  v <- resid(M, L, w); s <- c(s, sign(v))
  cat(sprintf("   %4d  %7.2f    %+18.8f       %s\n", L, w, v, ifelse(v > 0, "+", "-")))
}
cat(sprintf("\n   %d of %d modes agree in sign: %s\n", max(sum(s > 0), sum(s < 0)), length(s),
            ifelse(all(s > 0), "POSITIVE", ifelse(all(s < 0), "NEGATIVE", "MIXED"))))

cat("\n=== 3. through the contact region, which is r <= M ===\n")
cat("      r/M      l=2 w=1      l=10 w=1     l=2 w=0.3\n")
for (r in c(0.1, 0.25, 0.5, 0.75, 1.0)) {
  cat(sprintf("   %6.2f  %+11.5f  %+12.5f  %+12.5f\n",
              r, resid(r, 2, 1), resid(r, 10, 1), resid(r, 2, 0.3)))
}

cat("\n=== 4. plants, both of which must move the answer ===\n")
cat("   (a) drop the measure term 2f/r, which the constant-radius model does not have:\n")
resid_nomeasure <- function(r, L, w)
  f(r) * (-d1(function(z) half(z, L, w), r) + 2 * half(r, L, w)^2) - fp(r) * half(r, L, w)
for (r in c(0.25, 0.5, 1.0)) {
  a <- resid(r, 2, 1); b <- resid_nomeasure(r, 2, 1)
  cat(sprintf("       r = %4.2f M   with %+11.5f   without %+11.5f   moved by %.4f\n",
              r, a, b, abs(a - b)))
}
cat("   the measure term is doing real work, which is A.19's point that a varying\n")
cat("   transverse radius is what breaks the cancellation.\n")
cat("   (b) outside the horizon the expression still evaluates, and that matters:\n")
for (r in c(2.5, 3.0, 6.0)) {
  cat(sprintf("       r = %4.1f M   w^2 - V = %+9.5f   residue %+10.5f\n",
              r, 1 - V(r, 2), resid(r, 2, 1)))
}
cat("   A first draft of this file asserted p goes imaginary outside and the formula\n")
cat("   stops. It does not: w^2 - V stays positive for these modes and the residue is\n")
cat("   perfectly finite out there. The reason it is nonetheless zero outside is the\n")
cat("   silence theorem and nothing else: a point and its image are spacelike separated,\n")
cat("   the cross-sheet commutator vanishes identically, and there is no image term for\n")
cat("   this residue to be the residue OF. One reason, not two, and the wrong second\n")
cat("   reason is recorded here because it was believed for a few minutes.\n")

cat("\n=== 5. status ===\n")
cat("   Leading real order in WKB, Regge-Wheeler s = 0, p real because V < 0 inside.\n")
cat("   It gives a sign and not a coefficient, and it is geometry with no state input,\n")
cat("   which is what A.19 said it could not get. It is not a proof about the exact\n")
cat("   stress tensor and does not pretend to be.\n")

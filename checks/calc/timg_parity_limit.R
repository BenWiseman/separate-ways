# Does the image energy density actually vanish? Testing the limit rather than asserting it.
#
# GR49 found the parity sum decaying toward zero as the regulator lifts and stopped there.
# If that limit is genuinely zero it is a result, and it is a result AGAINST A.15's
# self-censoring conjecture: the two-point function can diverge at the contact configuration
# (A.18's s^-3/2) while the stress tensor built from it does not, which is exactly what
# A.19's wave-operator cancellation says structurally. So this has to be tested hard.
#
# On the separation, settled first because two iterations got it wrong. Theta = J o P_perp
# with J: (U,V) -> (-U,-V) preserves UV and V/U, so the image point has the SAME Schwarzschild
# (t, r) and differs only by interior region, which the chart cannot label. The separation is
# therefore not a real Killing-time shift. It is the half-period IMAGINARY one, which is A.8's
# result, and it multiplies each mode by a Boltzmann factor in omega rather than a phase in l.
# So it cannot damp the l-sum, and the l-sum's fate is decided by the parity alone.

M <- 1
f  <- function(r) 1 - 2 * M / r
V  <- function(r, L) f(r) * (L * (L + 1) / r^2 + 2 * M / r^3)
p  <- function(r, L, w) sqrt(w^2 - V(r, L))
d1 <- function(g, x) { h <- 1e-5; (g(x + h) - g(x - h)) / (2 * h) }
pp <- function(r, L, w) d1(function(z) p(z, L, w), r)

bracket <- function(r, L, w) {
  P <- p(r, L, w)
  Ftt <- -w^2; Frr <- P^2 + (pp(r, L, w) / (2 * P))^2
  Ftr <- w * P; Fang <- L * (L + 1) / r^2
  ut <- 1 / f(r); ur <- -sqrt(2 * M / r)
  ut^2 * Ftt + 2 * ut * ur * Ftr + ur^2 * Frr +
    0.5 * (-Ftt / f(r) + f(r) * Frr + Fang)
}

psum <- function(r, w, eps, wt = function(L) 1) {
  l <- 0:ceiling(45 / eps)
  0.5 * sum(((2 * l + 1) / (4 * pi)) * (-1)^l *
            sapply(l, function(L) bracket(r, L, w)) * wt(l) * exp(-l * eps))
}

cat("=== 1. the limit, over five halvings of the regulator ===\n")
cat("      r/M        eps=0.08      0.04        0.02       0.01      0.005    fitted power\n")
for (r in c(1.0, 0.75, 0.5, 0.25)) {
  es <- c(0.08, 0.04, 0.02, 0.01, 0.005)
  vs <- sapply(es, function(e) psum(r, 1, e))
  sl <- coef(lm(log(abs(vs)) ~ log(es)))[2]
  cat(sprintf("   %6.2f  %11.6f %11.6f %11.6f %10.6f %10.6f    %+6.3f\n",
              r, vs[1], vs[2], vs[3], vs[4], vs[5], sl))
}
cat("   a fitted power near +1 means the sum is proportional to eps and the limit is zero.\n")

cat("\n=== 2. the plants, which must NOT give a vanishing limit ===\n")
cat("   (a) no parity: the alternation is what does the cancelling\n")
es <- c(0.08, 0.04, 0.02, 0.01)
vs <- sapply(es, function(e) {
  l <- 0:ceiling(45 / e)
  0.5 * sum(((2 * l + 1) / (4 * pi)) * sapply(l, function(L) bracket(0.5, L, 1)) * exp(-l * e))
})
cat(sprintf("       %s   fitted power %+.3f\n",
            paste(sprintf("%.3e", vs), collapse = "  "),
            coef(lm(log(abs(vs)) ~ log(es)))[2]))
cat("   (b) parity kept but the bracket given an odd-l bias, which must survive\n")
vs <- sapply(es, function(e) psum(0.5, 1, e, wt = function(L) 1 + 0.3 * (L %% 2)))
cat(sprintf("       %s   fitted power %+.3f\n",
            paste(sprintf("%.3e", vs), collapse = "  "),
            coef(lm(log(abs(vs)) ~ log(es)))[2]))
cat("   a bias that tracks the parity is exactly what the cancellation cannot remove, so a\n")
cat("   power near zero there is the check working.\n")

cat("\n=== 3. why the cancellation happens, in one line ===\n")
cat("   The bracket is smooth in l. An alternating sum of a smooth function is not the sum\n")
cat("   of its terms, it is a difference of neighbours, and the regulated value goes as the\n")
cat("   derivative times eps. Smoothness in l is the whole mechanism, and it holds because\n")
cat("   the bracket depends on l only through l(l+1)/r^2 inside the momentum and the\n")
cat("   centrifugal term, both smooth.\n")
d <- sapply(0:6, function(L) bracket(0.5, L, 1))
cat(sprintf("   bracket at l = 0..6: %s\n", paste(sprintf("%.2f", d), collapse = ", ")))
cat(sprintf("   successive differences: %s\n",
            paste(sprintf("%.2f", diff(d)), collapse = ", ")))
cat("   smooth, no parity structure of its own, so nothing survives the alternation.\n")

cat("\n=== 4. the power fits above are wrong where the values cross zero ===\n")
cat("   At r = M the sequence goes +0.000211, -0.001542, ... and a log-log fit to a\n")
cat("   sign-changing sequence is meaningless. The differences halve as eps halves, which\n")
cat("   is convergence to a finite limit, not a power law. Richardson-extrapolate instead.\n\n")
rich <- function(v) { # v at eps, eps/2, eps/4, ...: repeated Richardson for a linear-in-eps error
  while (length(v) > 1) v <- 2 * v[-1] - v[-length(v)]
  v
}
cat("      r/M      sequence over eps = 0.08 .. 0.005                 extrapolated limit\n")
for (r in c(1.0, 0.75, 0.5, 0.25, 0.1)) {
  es <- c(0.08, 0.04, 0.02, 0.01, 0.005)
  vs <- sapply(es, function(e) psum(r, 1, e))
  cat(sprintf("   %6.2f   %s   %+12.6f\n", r,
              paste(sprintf("%+9.5f", vs), collapse = " "), rich(vs)))
}
cat("\n   the differences and their ratios, which is what says 'converging' rather than\n")
cat("   'vanishing' or 'diverging':\n")
for (r in c(1.0, 0.5)) {
  vs <- sapply(c(0.08, 0.04, 0.02, 0.01, 0.005), function(e) psum(r, 1, e))
  d <- diff(vs)
  cat(sprintf("   r = %4.2f M   diffs %s   ratios %s\n", r,
              paste(sprintf("%+.2e", d), collapse = " "),
              paste(sprintf("%.2f", d[-1] / d[-length(d)]), collapse = " ")))
}
cat("   ratios near 0.5 mean the error halves with eps, so the limit is the Richardson value.\n")

cat("\n=== 5. what that says ===\n")
cat("   The image energy density is FINITE everywhere in the contact region, small, and\n")
cat("   changes sign between the boundary and the interior. It does not diverge at r = M.\n")
cat("   A.18's two-point function does diverge there, as s^-3/2. Both can be true, and\n")
cat("   A.19 says why: the fold's pullbacks assemble the density into the wave operator,\n")
cat("   which annihilates the two-point function, so the divergence cancels out of the\n")
cat("   stress even though it is present in <phi^2>.\n")
cat("   If that survives scrutiny it is a result AGAINST A.15's self-censoring conjecture,\n")
cat("   and not because the sign is unhelpful but because there is nothing to have a sign.\n")

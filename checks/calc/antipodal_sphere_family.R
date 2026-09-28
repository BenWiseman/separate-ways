# The antipodal image two-point function on R x S^N, in closed form for every N, as a test
# of the caustic-order exponent rule the appendix now leans on.
#
# contact_world_function.py states: a Green function's generic light-cone singularity is
# steepened by half a power for every direction that refocuses on the connecting family. It
# rests on three points, one of them exact and two measured, and it is now load-bearing. This
# tries to break it on a family where the answer is exact at every member.
#
# THE FAMILY. On R x S^N of radius a, the conformally coupled massless scalar has
# xi = (N-1)/(4N) in D = N+1 dimensions and R = N(N-1)/a^2, so
#     omega_n^2 = [n(n + N - 1) + (N-1)^2/4]/a^2 = [n + (N-1)/2]^2/a^2,
# equally spaced, which is the whole reason these sums close. Write lambda = (N-1)/2. The
# Gegenbauer addition theorem gives the mode sum at angle gamma in terms of C_n^lambda, and
# at the ANTIPODE C_n^lambda(-1)/C_n^lambda(1) = (-1)^n, so the parity arrives on its own
# exactly as it does at N = 3.

Vol <- function(N) 2 * pi^((N + 1) / 2) / gamma((N + 1) / 2)
lam <- function(N) (N - 1) / 2
# in log space, or gamma() overflows well before the sums have converged
deg <- function(n, N) {
  L <- lam(N)
  exp(log(2 * (n + L)) + lgamma(n + 2 * L) - lgamma(2 * L + 1) - lgamma(n + 1))
}

cat("=== 1. the degeneracies are the familiar ones ===\n")
cat("      N   Vol(S^N)      d_0  d_1  d_2  d_3     known\n")
known <- list("2" = "2n+1", "3" = "(n+1)^2", "4" = "(n+1)(n+2)(2n+3)/6")
for (N in 2:5) {
  d <- sapply(0:3, deg, N = N)
  cat(sprintf("   %4d  %10.6f   %4.0f %4.0f %4.0f %4.0f     %s\n",
              N, Vol(N), d[1], d[2], d[3], d[4],
              if (!is.null(known[[as.character(N)]])) known[[as.character(N)]] else ""))
}
stopifnot(all(abs(sapply(0:5, deg, N = 2) - (2 * (0:5) + 1)) < 1e-9))
stopifnot(all(abs(sapply(0:5, deg, N = 3) - ((0:5) + 1)^2) < 1e-9))

cat("\n=== 2. the sum closes, for every N ===\n")
cat("   (1/2 omega_n) d_n cancels the (n + lambda), leaving a binomial series:\n")
cat("     G_img = a^{1-N} / ( 2 lambda 2^{2 lambda} Vol(S^N) cos^{2 lambda}(eta/2) ).\n")
G_cf <- function(eta, N, a = 1) {
  L <- lam(N)
  a^(1 - N) / (2 * L * 2^(2 * L) * Vol(N) * cos(eta / 2)^(2 * L))
}
G_abel <- function(eta, N, a, K, damp) {
  L <- lam(N); n <- 0:K
  Re(sum((a / (2 * (n + L))) * deg(n, N) * (-1)^n / (a^N * Vol(N)) *
         exp(-1i * (n + L) * eta) * exp(-damp * n)))
}
# The damping biases the answer at O(damp), and more strongly the faster d_n grows, so take
# two dampings and extrapolate rather than quoting one.
# Three dampings and a quadratic Richardson: the bias is O(damp) only to leading order and
# the higher terms bite once d_n grows like n^4. K is set so the truncation is below e^-100.
G_sum <- function(eta, N, a = 1, K = 20000, damp = 0.02) {
  g1 <- G_abel(eta, N, a, K, damp)
  g2 <- G_abel(eta, N, a, K, damp / 2)
  g3 <- G_abel(eta, N, a, K, damp / 4)
  (8 * g3 - 6 * g2 + g1) / 3
}
cat("\n      N     eta      Abel sum        closed form     difference\n")
worst <- 0
for (N in 2:5) for (e in c(0.6, 1.4, 2.2)) {
  s <- G_sum(e, N); c0 <- G_cf(e, N); worst <- max(worst, abs(s - c0) / abs(c0))
  cat(sprintf("   %4d  %6.2f  %14.8f  %14.8f   %.2e\n", N, e, s, c0, abs(s - c0)))
}
cat(sprintf("   worst RELATIVE difference over the twelve points: %.1e\n", worst))
stopifnot(worst < 1e-3)

cat("\n=== 3. N = 3 must reproduce the appendix's own closed form ===\n")
for (e in c(1.0, 2.0, 3.0)) {
  cat(sprintf("   eta = %.1f:  general formula %14.9f   1/(16 pi^2 cos^2(eta/2)) = %14.9f\n",
              e, G_cf(e, 3), 1 / (16 * pi^2 * cos(e / 2)^2)))
  stopifnot(abs(G_cf(e, 3) - 1 / (16 * pi^2 * cos(e / 2)^2)) < 1e-12)
}
cat("   So the N = 3 member is the Einstein static universe and the family contains it.\n")

cat("\n=== 4. the exponent, and the rule it is supposed to obey ===\n")
cat("   Near eta = pi, cos(eta/2) ~ delta/2, so G_img ~ delta^{-2 lambda} = delta^{-(N-1)}.\n")
cat("   The connecting geodesics from a point to its antipode on S^N form an S^{N-1}, so the\n")
cat("   caustic order is n = N - 1 and the spacetime dimension is D = N + 1. The rule as\n")
cat("   stated in contact_world_function.py is -(1 + n/2), which is a FOUR-dimensional\n")
cat("   specialisation. The general statement has to reduce to -(D-2)/2 with no caustic:\n")
cat("        exponent = -(D - 2 + n)/2.\n\n")
cat("      N   D   n     measured exponent   -(D-2+n)/2    -(1 + n/2)\n")
for (N in 2:6) {
  D <- N + 1; n <- N - 1
  d1 <- 1e-3; d2 <- 1e-4
  p <- log(G_cf(pi - d2, N) / G_cf(pi - d1, N)) / log(d2 / d1)
  cat(sprintf("   %4d %3d %3d   %17.6f   %10.1f   %11.1f\n", N, D, n, p, -(D - 2 + n)/2, -(1 + n/2)))
  stopifnot(abs(p - (-(D - 2 + n) / 2)) < 1e-6)
}
cat("   The general rule holds at every member; the four-dimensional shorthand does not,\n")
cat("   and agrees only where it should, at D = 4.\n")

cat("\n=== 5. the rule against everything else that is known ===\n")
cat("      case                                   D   n    rule      known        source\n")
rows <- list(
  list("flat space, no caustic",                   4, 0, -1.0,     "measured -1.0095", "A.18 control"),
  list("A.18 image pair",                          4, 1, -1.5,     "measured -1.49874", "A.18"),
  list("Einstein static universe",                 4, 2, -2.0,     "exact -2",          "section 3 above"),
  list("R x S^4",                                  5, 3, -3.0,     "exact -3",          "section 4 above"),
  list("R x S^5",                                  6, 4, -4.0,     "exact -4",          "section 4 above"))
for (r in rows) {
  rule <- -(r[[2]] - 2 + r[[3]]) / 2
  cat(sprintf("      %-38s %2d  %2d  %7.3f   %-18s %s\n", r[[1]], r[[2]], r[[3]], rule, r[[5]], r[[6]]))
  stopifnot(abs(rule - r[[4]]) < 1e-9)
}
cat("   Five points spanning three spacetime dimensions and five caustic orders, two of\n")
cat("   them measured and three exact, and none of them fitted to the rule.\n")

cat("\n=== 6. what it says about a black hole ===\n")
cat("   contact_conjugacy_general.R gives n = D - 3 for a spherically symmetric hole, one\n")
cat("   fewer than this family because the d_t direction spreads instead of focusing. So\n")
cat("        exponent = -(D - 2 + D - 3)/2 = -(2D - 5)/2,\n")
cat("   which is -3/2 in four dimensions and the stress two powers steeper at -7/2.\n")
for (D in 4:7) cat(sprintf("      D = %d:  n = %d, image term delta^{%+.1f}, stress delta^{%+.1f}\n",
                           D, D - 3, -(2 * D - 5) / 2, -(2 * D - 5) / 2 - 2))

cat("\n=== 7. the sign, which comes free ===\n")
cat("   cos(eta/2) > 0 for eta in (0, pi), so G_img > 0 at every N whatever the parity of\n")
cat("   2 lambda. The image two-point function approaches contact from ABOVE in every\n")
cat("   dimension, which is the sign A.18 measured and the sign its 502 carries.\n")
for (N in 2:5) cat(sprintf("      N = %d:  G_img at eta = pi - 0.01 is %+.4e\n", N, G_cf(pi - 0.01, N)))
cat("   That is the two-point function's sign and not the stress tensor's, which\n")
cat("   image_stress_conformal.R shows can be zero while this is positive.\n")

cat("\n=== 8. the plant ===\n")
cat("   (a) remove the antipodal parity and the divergence must leave the antipode.\n")
G_nopar <- function(eta, N, K = 20000, damp = 0.02) {
  L <- lam(N); n <- 0:K
  f <- function(dp) Re(sum((1 / (2 * (n + L))) * deg(n, N) / Vol(N) *
                           exp(-1i * (n + L) * eta) * exp(-dp * n)))
  (8 * f(damp / 4) - 6 * f(damp / 2) + f(damp)) / 3
}
for (N in c(2, 3)) {
  cat(sprintf("      N = %d:  with parity at eta = pi - 0.05: %12.3f   without: %12.3f\n",
              N, G_sum(pi - 0.05, N), G_nopar(pi - 0.05, N)))
  cat(sprintf("               with parity at eta = 0.05:      %12.3f   without: %12.3f\n",
              G_sum(0.05, N), G_nopar(0.05, N)))
}
cat("   The parity moves the divergence from coincidence to the antipode, at every N.\n")
cat("   (b) the exponent estimator must report the right thing on a known control.\n")
for (p0 in c(-1, -1.5, -3)) {
  f <- function(d) d^p0
  est <- log(f(1e-4) / f(1e-3)) / log(1e-4 / 1e-3)
  cat(sprintf("      planted exponent %+5.2f  ->  estimator returns %+8.5f\n", p0, est))
  stopifnot(abs(est - p0) < 1e-9)
}

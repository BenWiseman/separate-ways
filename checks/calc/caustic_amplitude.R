# The caustic AMPLITUDE, which is the last thing the exponent rule does not give.
#
# antipodal_sphere_family.R fixed the exponent at -(D-2+n)/2 on five points. The coefficient
# is what decides the magnitude of A.15's self-censoring stress, and the exponent argument
# says nothing about it. This gets at it through the heat kernel, where a caustic shows up
# not as a divergence but as a replacement:
#
#     G = int_0^inf ds (4 pi s)^{-D/2} Delta^{1/2} e^{-sigma/2s},
#
# and at an n-fold degenerate connecting family Delta^{1/2} is infinite. What replaces it is
# finite and carries s:
#
#     Delta^{1/2}  ->  c_n s^{-n/2},
#
# which is exactly what turns sigma^{-(D-2)/2} into sigma^{-(D-2+n)/2}. So the exponent rule
# and the amplitude are the same statement, and c_n is the number wanted.

a <- 1
Vol  <- function(N) 2 * pi^((N + 1) / 2) / gamma((N + 1) / 2)
deg  <- function(n, N) exp(log(2 * (n + (N - 1) / 2)) + lgamma(n + N - 1) -
                           lgamma(N) - lgamma(n + 1))

cat("=== 1. the prediction, read off the exact Green functions ===\n")
cat("   Feeding Delta^{1/2} = c_n s^{-n/2} through the proper-time integral gives\n")
cat("     G = c_n (4 pi)^{-D/2} Gamma((D+n)/2 - 1) (2/sigma)^{(D+n)/2 - 1},\n")
cat("   and matching that against G_img = a^{1-N}/(2 lambda Vol(S^N) delta^{N-1}) with\n")
cat("   sigma = pi a^2 delta, D = N+1 and n = N-1 leaves\n")
cat("     c_n = L^n pi^{n/2} 4 pi / ( n! Vol(S^{n+1}) ),   L = pi a the geodesic length.\n\n")
c_pred <- function(N) { n <- N - 1; L <- pi * a; L^n * pi^(n/2) * 4 * pi / (factorial(n) * Vol(N)) }
cat("      N   n   c_n predicted        in closed form\n")
forms <- c("pi^{3/2} a", "pi^2 a^2", "pi^{7/2} a^3 / 4", "pi^4 a^4 / 6")
closed <- c(pi^1.5, pi^2, pi^3.5 / 4, pi^4 / 6)
for (N in 2:5) {
  cat(sprintf("   %4d %3d   %16.9f     %-18s %16.9f\n",
              N, N - 1, c_pred(N), forms[N - 1], closed[N - 1]))
  stopifnot(abs(c_pred(N) - closed[N - 1]) < 1e-9)
}

cat("\n=== 2. the heat kernel at the antipode, independently ===\n")
cat("   K(pi; s) = (1/(a^N Vol(S^N))) sum_n d_n (-1)^n exp(-n(n+N-1) s/a^2), the parity again\n")
cat("   arriving from the addition theorem rather than being inserted.\n")
Kanti <- function(s, N, K = 400) {
  n <- 0:K
  sum(deg(n, N) * (-1)^n * exp(-n * (n + N - 1) * s / a^2)) / (a^N * Vol(N))
}
cat("   S^3 has an exact closed form to check the sum against:\n")
cat("     K = (4 pi s)^{-3/2} e^{s} 2 e^{-pi^2/4s} (pi^2/2s - 1)\n")
K3_cf <- function(s) (4 * pi * s)^(-1.5) * exp(s) * 2 * exp(-pi^2 / (4 * s)) * (pi^2 / (2 * s) - 1)
cat("        s        mode sum         closed form      rel. diff\n")
for (s in c(0.5, 0.3, 0.2, 0.15)) {
  m <- Kanti(s, 3); c0 <- K3_cf(s)
  cat(sprintf("   %6.3f  %15.8e  %15.8e   %.2e\n", s, m, c0, abs(m / c0 - 1)))
  stopifnot(abs(m / c0 - 1) < 1e-6)
}
cat("   So the mode sum is trustworthy where it is used.\n")

cat("\n=== 3. extracting c_n, with the subleading terms fitted rather than ignored ===\n")
cat("   A(s) = K(pi;s) (4 pi s)^{N/2} e^{pi^2 a^2/4s} s^{n/2} tends to c_n as s -> 0, with\n")
cat("   an expansion in s on the way, so fit and extrapolate to zero. Orthogonal polynomials,\n")
cat("   degree five, over s in [0.10, 0.45]: below 0.10 the alternating sum loses digits to\n")
cat("   cancellation and the fit degrades, which section 6 checks rather than assumes.\n\n")
cat("      N   n   c_n from the heat kernel   predicted      rel. diff   fit residual\n")
for (N in 2:5) {
  n <- N - 1
  ss <- seq(0.10, 0.45, length.out = 20)
  A  <- sapply(ss, function(s) Kanti(s, N) * (4 * pi * s)^(N / 2) * exp(pi^2 * a^2 / (4 * s)) * s^(n / 2))
  f  <- lm(A ~ poly(ss, 5))
  c0 <- unname(predict(f, newdata = data.frame(ss = 0))); res <- max(abs(residuals(f)))
  cat(sprintf("   %4d %3d   %22.8f   %12.8f   %.2e   %.1e\n",
              N, n, c0, c_pred(N), abs(c0 / c_pred(N) - 1), res))
  stopifnot(abs(c0 / c_pred(N) - 1) < 1e-5)
}
cat("\n   Two independent routes, a Lorentzian Green function and a Riemannian heat kernel,\n")
cat("   land on the same four numbers. The caustic amplitude is c_n s^{-n/2} and c_n is\n")
cat("   L^n pi^{n/2} 4 pi / (n! Vol(S^{n+1})) on this family.\n")

cat("\n=== 4. the amplitude does not know about flat directions ===\n")
cat("   The heat kernel on a metric product factorises, K = K_1 K_2, and a flat factor\n")
cat("   contributes (4 pi s)^{-k/2} with Delta = 1. So Delta^{1/2} for M x R^k equals\n")
cat("   Delta^{1/2} for M, and c_n is a property of the degenerate family alone.\n")
for (k in 0:2) {
  s <- 0.2; N <- 2; D <- N + k
  Kprod <- Kanti(s, N) * (4 * pi * s)^(-k / 2)
  A <- Kprod * (4 * pi * s)^(D / 2) * exp(pi^2 / (4 * s)) * s^(1 / 2)
  cat(sprintf("      S^2 x R^%d (D = %d):  A(s = 0.2) = %.8f\n", k, D, A))
}
cat("   Identical, as it must be. A.18's geometry is R^2 x S^2, which is the k = 2 row, so\n")
cat("   ITS caustic amplitude is S^2's: n = 1 in D = 4, the case a black hole is in.\n")
cat(sprintf("      c_1 = pi^{3/2} a = %.8f a, with a the transverse radius.\n", pi^1.5))

cat("\n=== 5. what this does NOT give, stated plainly ===\n")
cat("   c_1 = pi^{3/2} a on S^2, and L = pi a there, so c_1 = sqrt(pi) L. Which of those is\n")
cat("   the right way to write it cannot be decided from one example: on a round sphere\n")
cat("   every length is proportional to a, so L, the maximum of the Jacobi field, and the\n")
cat("   integral of it are all the same number up to a constant. A black hole's contact\n")
cat("   geodesic has a DIFFERENT Jacobi profile, M(1 + sin phi) sin phi against a sin(u/a),\n")
cat("   so the three candidates separate there and give different answers:\n")
Lsch <- 3 * pi / 2 + 4
Jmax <- max(sapply(seq(0, pi, length.out = 20001), function(p) (1 + sin(p)) * sin(p)))
Jint <- integrate(function(p) (1 + sin(p)) * sin(p) * (1 + sin(p))^2, 0, pi)$value  # d(lambda) = r^2 dphi
cat(sprintf("      affine length L                              = %.6f M\n", Lsch))
cat(sprintf("      max of the focusing Jacobi field             = %.6f M\n", Jmax))
cat(sprintf("      int J dlambda over the curve                 = %.6f M^2\n", Jint))
cat(sprintf("      sqrt(pi) L would give                        = %.6f M\n", sqrt(pi) * Lsch))
cat("   A second n = 1 example with a different profile would separate them. That example now\n")
cat("   exists and this section is superseded by it: caustic_amplitude_profile.R builds a\n")
cat("   non-symmetric surface of revolution, which moves the Jacobi profile while holding the\n")
cat("   length and the family volume fixed, and the heat kernel between the poles is unmoved. So\n")
cat("   the amplitude does not see the profile at all, which leaves the affine length: c_1 =\n")
cat("   sqrt(pi) L, the last line of the table above, and at Schwarzschild 15.442307 M. Do not\n")
cat("   read the sentence this section used to end with, that the amplitude is undetermined. It\n")
cat("   was true when written and it is not now.\n")

cat("\n=== 6. the plants ===\n")
cat("   (a) the extraction must fail if the power is wrong.\n")
for (nwrong in 0:4) {
  N <- 3; ss <- seq(0.10, 0.45, length.out = 14)
  A <- sapply(ss, function(s) Kanti(s, N) * (4 * pi * s)^(N / 2) * exp(pi^2 / (4 * s)) * s^(nwrong / 2))
  sp <- diff(range(A)) / abs(mean(A))
  cat(sprintf("      S^3, dividing by s^{-%d/2}:  A varies by %8.1f%% over the window%s\n",
              nwrong, 100 * sp, if (nwrong == 2) "   <- the minimum" else ""))
  spread <- c(if (exists("spread")) spread, sp)
}
stopifnot(which.min(spread) == 3)
cat("   The variation is a MINIMUM at n = 2 and rises on both sides, so the power is\n")
cat("   measured. A test that only fell monotonically would prove nothing.\n")
cat("   (a2) and the window's lower edge is where it is for a reason: below 0.10 the\n")
cat("        alternating sum cancels away its own digits.\n")
for (lo in c(0.10, 0.08, 0.06)) {
  N <- 3; ss <- seq(lo, 0.45, length.out = 20)
  A <- sapply(ss, function(s) Kanti(s, N) * (4 * pi * s)^(N / 2) * exp(pi^2 / (4 * s)) * s^(1))
  c0 <- unname(predict(lm(A ~ poly(ss, 5)), newdata = data.frame(ss = 0)))
  cat(sprintf("        window from %.2f:  c_2 = %14.6f   error %.1e\n",
              lo, c0, abs(c0 / c_pred(3) - 1)))
}
cat("   (b) remove the antipodal parity and the heat kernel at the antipode must change.\n")
Kno <- function(s, N, K = 400) { n <- 0:K; sum(deg(n, N) * exp(-n * (n + N - 1) * s / a^2)) / (a^N * Vol(N)) }
for (N in c(2, 3)) cat(sprintf("      N = %d, s = 0.2:  with parity %.4e   without %.4e\n",
                               N, Kanti(0.2, N), Kno(0.2, N)))
cat("   Without it the sum is the coincidence-limit kernel and has no antipodal suppression.\n")

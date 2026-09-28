# The S^4 antipodal propagator, and the four partial sums C.1 quotes.
#
# Appendix C.1 states that the conformal propagator at the antipode, G = 1/(16 pi^2 R^2), is
# "obtained also as an Abel-summed spectral series with partial sums 1.08033, 1.00755, 1.00075,
# 1.000075". Those four numbers had no script behind them: an audit of every high-precision
# number in the cosmology paper against every script in both release trees found them reproduced
# by nothing. They are correct, and this is the file that says so.
#
# THE SERIES. On S^4 of radius R the conformally coupled scalar has xi = 1/6 and Ricci scalar
# 12/R^2, so its eigenvalues are [n(n+3) + 2]/R^2 = (n+1)(n+2)/R^2 with degeneracy
# d_n = (n+1)(n+2)(2n+3)/6. The addition theorem at the antipode contributes (-1)^n d_n / Vol,
# and d_n over the eigenvalue collapses to R^2 (2n+3)/6, so with Vol(S^4) = 8 pi^2 R^4/3
#
#   G(x, Ax) = (1/16 pi^2 R^2) sum_n (-1)^n (2n+3),
#
# a divergent alternating series whose Abel sum is 1. The quoted figures are that Abel sum at
# damping x = 0.9, 0.99, 0.999 and 0.9999.

R <- 1
deg <- function(n) (n+1)*(n+2)*(2*n+3)/6
eig <- function(n) (n+1)*(n+2)/R^2
VolS4 <- 8*pi^2*R^4/3

cat("=== 1. the spectrum, and the collapse of degeneracy over eigenvalue ===\n")
cat("      n    eigenvalue (n+1)(n+2)   degeneracy   ratio d_n/eig    R^2(2n+3)/6\n")
for (n in 0:5)
  cat(sprintf("   %4d   %19.0f   %10.0f   %13.6f   %13.6f\n",
              n, eig(n), deg(n), deg(n)/eig(n), R^2*(2*n+3)/6))
stopifnot(max(abs(sapply(0:60, function(n) deg(n)/eig(n) - R^2*(2*n+3)/6))) < 1e-9)
cat("   Exact at every n to 60, so the series really is sum (-1)^n (2n+3) up to the prefactor.\n")

cat("\n=== 2. the Abel sum, in closed form and by direct summation ===\n")
cat("   sum_n (-1)^n (2n+3) x^n = 3/(1+x) - 2x/(1+x)^2 = (3+x)/(1+x)^2, which is 1 at x = 1.\n\n")
closed <- function(x) (3 + x)/(1 + x)^2
direct <- function(x, N = 400000) sum((-1)^(0:N) * (2*(0:N) + 3) * x^(0:N))
cat("        x          closed form      direct sum        difference\n")
quoted <- c(1.08033, 1.00755, 1.00075, 1.000075)
xs <- c(0.9, 0.99, 0.999, 0.9999)
for (i in seq_along(xs)) {
  cf <- closed(xs[i]); ds <- direct(xs[i])
  cat(sprintf("   %8.4f   %15.8f   %14.8f   %14.2e\n", xs[i], cf, ds, abs(cf - ds)))
  stopifnot(abs(cf - ds) < 1e-6)
}
cat(sprintf("\n   limit at x = 1: %.10f\n", closed(1)))
stopifnot(abs(closed(1) - 1) < 1e-14)

cat("\n=== 3. against the four figures C.1 quotes ===\n")
# the appendix quotes each to its own number of decimals, so the bar is half a unit in the last
# place quoted and not a fixed tolerance; a flat 5e-7 fails on 1.08033 for rounding alone
dp <- function(q) { t <- format(q, scientific = FALSE); nchar(sub(".*\\.", "", t)) }
cat("        x          computed        quoted     decimals    half a unit    agreement\n")
for (i in seq_along(xs)) {
  cf <- closed(xs[i]); d <- dp(quoted[i]); tol <- 0.5 * 10^(-d)
  cat(sprintf("   %8.4f   %14.7f   %11.6f   %8d   %12.1e   %11.1e\n",
              xs[i], cf, quoted[i], d, tol, abs(cf - quoted[i])))
  stopifnot(abs(cf - quoted[i]) <= tol)
}
cat("   All four reproduce, so the appendix's series is arithmetic and not recollection.\n")

cat("\n=== 4. and the closed form it is checked against ===\n")
G_closed <- 1/(16*pi^2*R^2)
G_series <- (1/(16*pi^2*R^2)) * closed(1)
cat(sprintf("   G(x, Ax) closed form      = %.12f\n", G_closed))
cat(sprintf("   from the series at x -> 1 = %.12f\n", G_series))
cat(sprintf("   difference                = %.2e\n", abs(G_closed - G_series)))
cat("   which is what C.1 means by the agreement checking an implementation rather than proving\n")
cat("   anything independent: both sides are the same fixed identity, reached two ways.\n")
stopifnot(abs(G_closed - G_series) < 1e-15)

cat("\n=== 5. the plant: the antipodal parity is what makes the sum finite ===\n")
cat("   Drop the (-1)^n, which is the same as asking for the propagator at coincidence rather\n")
cat("   than at the antipode, and the Abel sum must blow up instead of approaching one.\n\n")
nopar <- function(x) sum((2*(0:200000) + 3) * x^(0:200000))
cat("        x          with parity     without parity\n")
for (x in c(0.9, 0.99, 0.999))
  cat(sprintf("   %8.4f   %14.6f   %16.2f\n", x, closed(x), nopar(x)))
cat("   Without the parity it diverges as x approaches one, which is the coincidence singularity\n")
cat("   of a two-point function sitting where it belongs. The parity moves it to the antipode.\n")
stopifnot(nopar(0.999) > 1e5, closed(0.999) < 1.01)

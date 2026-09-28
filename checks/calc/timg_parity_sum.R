# The hole in GR48: a positive bracket is not a positive sum.
#
# A.19's density is  rho = (1/2) sum_l c_l (-1)^l [ bracket_l ],  c_l = (2l+1)/4pi, and the
# (-1)^l is the antipodal parity, the same factor A.18 shows IS the transverse map. GR48
# established that the bracket is positive at every mode tested and concluded the density
# is positive. That does not follow. With the parity alternating, a bracket that is positive
# and growing in l gives a sum whose sign is decided by the alternation, not by the bracket.
#
# The terms grow: GR48 has +6.4, +24.8, +254 at l = 0, 2, 10, and c_l ~ l, so the summand
# goes as l^3 with alternating sign. That is A.18's situation exactly, and it needs A.18's
# regulator: the Wightman function's own i-epsilon, which damps the tower by exp(-l eps),
# with eps taken small and the answer checked for being free of it.

M <- 1
f  <- function(r) 1 - 2 * M / r
V  <- function(r, L) f(r) * (L * (L + 1) / r^2 + 2 * M / r^3)
p  <- function(r, L, w) sqrt(w^2 - V(r, L))
d1 <- function(g, x) { h <- 1e-5; (g(x + h) - g(x - h)) / (2 * h) }
pp <- function(r, L, w) d1(function(z) p(z, L, w), r)

F_tt  <- function(r, L, w) -w^2
F_rr  <- function(r, L, w) { P <- p(r, L, w); P^2 + (pp(r, L, w) / (2 * P))^2 }
F_tr  <- function(r, L, w) w * p(r, L, w)
F_ang <- function(r, L, w) L * (L + 1) / r^2

bracket <- function(r, L, w) {                 # the infalling density, per mode
  ut <- 1 / f(r); ur <- -sqrt(2 * M / r)
  sq <- ut^2 * F_tt(r, L, w) + 2 * ut * ur * F_tr(r, L, w) + ur^2 * F_rr(r, L, w)
  tr <- -F_tt(r, L, w) / f(r) + f(r) * F_rr(r, L, w) + F_ang(r, L, w)
  sq + 0.5 * tr
}

parity_sum <- function(r, w, eps, lmax = NULL) {
  if (is.null(lmax)) lmax <- ceiling(40 / eps)
  l <- 0:lmax
  c_l <- (2 * l + 1) / (4 * pi)
  b   <- sapply(l, function(L) bracket(r, L, w))
  0.5 * sum(c_l * (-1)^l * b * exp(-l * eps))
}

cat("=== 1. the bracket is positive and grows, which is what makes the sum non-obvious ===\n")
cat("      l        bracket at r = M, omega = 1\n")
for (L in c(0, 1, 2, 3, 10, 30, 100)) cat(sprintf("   %5d   %+18.4f\n", L, bracket(1, L, 1)))
cat("   growing roughly as l^2, and c_l ~ l, so the summand goes as l^3, alternating.\n")

cat("\n=== 2. the regulated alternating sum, and whether it settles ===\n")
cat("      eps        sum at r = M, omega = 1\n")
vals <- c()
for (e in c(0.2, 0.1, 0.05, 0.025, 0.0125)) {
  v <- parity_sum(1, 1, e); vals <- c(vals, v)
  cat(sprintf("   %8.4f   %+20.6f\n", e, v))
}
cat("   ratio of successive values (1 means converged, 2 means it is doubling):\n")
cat(sprintf("     %s\n", paste(sprintf("%.3f", vals[-1] / vals[-length(vals)]), collapse = "  ")))

cat("\n=== 3. the same at two more radii inside the contact region ===\n")
for (r in c(0.5, 0.25)) {
  v <- sapply(c(0.1, 0.05, 0.025), function(e) parity_sum(r, 1, e))
  cat(sprintf("   r = %5.2f M   eps 0.1, 0.05, 0.025:  %+12.4f  %+12.4f  %+12.4f\n",
              r, v[1], v[2], v[3]))
}

cat("\n=== 4. the control: drop the parity and the sum must behave completely differently ===\n")
nosum <- function(r, w, eps, lmax = ceiling(40 / eps)) {
  l <- 0:lmax
  0.5 * sum(((2 * l + 1) / (4 * pi)) * sapply(l, function(L) bracket(r, L, w)) * exp(-l * eps))
}
for (e in c(0.1, 0.05, 0.025)) {
  cat(sprintf("   eps %6.4f   with parity %+14.4f    without %+18.4f\n",
              e, parity_sum(1, 1, e), nosum(1, 1, e)))
}
cat("   without the parity the sum blows up like 1/eps^4, with it the alternation cancels\n")
cat("   most of the tower. That difference is the whole content of the question.\n")

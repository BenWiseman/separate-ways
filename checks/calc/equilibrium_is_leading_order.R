#!/usr/bin/env Rscript
# Equilibrium is not a separate assumption. It is the leading term of any Hadamard state.
#
# Section 3 says the squeeze at a horizon is fixed "subject to equilibrium being assumed rather
# than proved", and lists equilibrium as one of two costs. That entry can be paid off with
# something the construction already carries, because the image stress is computed from the
# Hadamard parametrix and so every state in play is Hadamard by assumption already.
#
# The Hadamard condition says a state's two-point function splits as W = W_sing + W_reg with
# W_sing fixed by the geometry alone and W_reg smooth. The modular temperature at a horizon comes
# from the short-distance structure, which is the part the condition makes state-independent. So
# the question is not whether a given state is in equilibrium but how much of the balance the
# state-dependent part can reach, and the answer is a power of the patch size.
#
# Measured here on flat four-dimensional space with a massless scalar, for two families of
# Hadamard states that are nowhere near the vacuum:
#
#   thermal at temperature T:  W_T - W_0 = (1/4 pi^2 beta^2)[pi coth(pi x)/x - 1/x^2], x = r/beta
#   coherent on a background:  W_c - W_0 = phi_cl(x) phi_cl(y)
#
# In both the ratio W_reg/W_sing vanishes as r^2 with a coefficient set by the state's own scale,
# so a patch of size r sees the vacuum to relative accuracy (r times that scale)^2. The Clausius
# balance is already the leading order in r, which puts the correction one order below the
# equation it returns. That is the same shape as the universality result: an isometry through
# second order, failing at third, with the failure not reaching the balance.
#
# The plant is the state that is NOT Hadamard, with the singular coefficient altered by delta. Then
# the ratio tends to delta instead of zero, at exponent zero, and the argument fails as it should.

TOL <- 1e-9
fail <- 0
note <- function(ok, what) {
  if (!ok) { cat("   *** FAILED:", what, "\n"); fail <<- fail + 1 }
  invisible(ok)
}

Wsing <- function(r) 1 / (4 * pi^2 * r^2)              # equal-time vacuum, 4D massless

# W_T - W_0 at equal times, from the image sum over Matsubara copies, in closed form. Written as
# [pi x coth(pi x) - 1] / x^2 rather than pi coth(pi x)/x - 1/x^2: the second form subtracts 1/x^2
# from itself and at x = 1e-6 that leaves four significant figures out of sixteen. Below pi x = 1/2
# the series is used instead, which is where the cancellation bites.
dW_thermal <- function(r, beta) {
  x <- r / beta
  u <- pi * x
  big <- (u * cosh(u) - sinh(u)) / sinh(u)           # = u coth(u) - 1, stable for u not small
  ser <- u^2 / 3 - u^4 / 45 + 2 * u^6 / 945 - u^8 / 4725
  num <- ifelse(u < 0.5, ser, big)
  (1 / (4 * pi^2 * beta^2)) * num / x^2
}

cat("=== 1. the closed form against the image sum it came from ===\n")
cat("   W_T = (1/4 pi^2) sum_n 1/[r^2 - (t - i n beta)^2]; at t = 0 the n =/= 0 terms are\n")
cat("   (1/2 pi^2) sum_{n>=1} 1/(r^2 + n^2 beta^2). Summed to 200000 terms with the tail\n")
cat("   integrated, against the closed form:\n")
cat("      r/beta        image sum          closed form        relative\n")
for (xr in c(0.05, 0.2, 1, 3)) {
  beta <- 1.7; r <- xr * beta
  # The tail is done in closed form, integral du/(r^2 + u^2 beta^2) = atan(u beta / r)/(r beta).
  # R's integrate() returns zero on this one, because with a lower limit of 2e5 the integrand is
  # 4e-12 there and the adaptive rule sees nothing: the whole 1.7e-6 tail went missing and the
  # comparison was off by exactly that.
  N <- 200000; n <- 1:N
  tail <- (pi / 2 - atan((N + 0.5) * beta / r)) / (r * beta)
  s <- sum(1 / (r^2 + n^2 * beta^2)) + tail
  got <- s / (2 * pi^2); want <- dW_thermal(r, beta)
  cat(sprintf("   %9.3f %18.10f %18.10f %14.2e\n", xr, got, want, abs(got / want - 1)))
  note(abs(got / want - 1) < 1e-6, "image sum reproduces the closed form")
}
cat("   And the coincidence limit is the textbook <phi^2>_T - <phi^2>_0 = T^2/12:\n")
for (beta in c(0.4, 1, 3.3)) {
  lim <- dW_thermal(1e-6, beta)
  cat(sprintf("      beta = %6.3f   limit = %14.10f   T^2/12 = %14.10f\n",
              beta, lim, 1 / (12 * beta^2)))
  note(abs(lim / (1 / (12 * beta^2)) - 1) < 1e-8, "coincidence limit is T^2/12")
}

cat("\n=== 2. the state-dependent part is quadratically small in the patch ===\n")
cat("   The ratio W_reg/W_sing over five decades of r/beta, with the exponent fitted on the\n")
cat("   smallest three decades and the coefficient compared with pi^2/3:\n")
ratio_thermal <- function(r, beta) dW_thermal(r, beta) / Wsing(r)
beta <- 1.0
rs <- 10^seq(-6, -1, length.out = 26)
rat <- ratio_thermal(rs, beta)
cat("      r/beta           W_reg/W_sing        that over r^2\n")
for (i in c(1, 6, 11, 16, 21, 26)) {
  cat(sprintf("   %12.3e %20.10e %20.10f\n", rs[i], rat[i], rat[i] / rs[i]^2))
}
use <- rs < 1e-3
ex <- coef(lm(log(rat[use]) ~ log(rs[use])))[2]
cf <- rat[1] / rs[1]^2
cat(sprintf("   fitted exponent %.6f against 2, coefficient %.7f against pi^2/3 = %.7f\n",
            ex, cf, pi^2 / 3))
note(abs(ex - 2) < 1e-5, "thermal exponent is 2")
note(abs(cf / (pi^2 / 3) - 1) < 1e-8, "thermal coefficient is pi^2/3")

cat("\n=== 3. the same for a coherent state, where the regular part is a classical field ===\n")
cat("   W_c - W_0 = phi_cl(x) phi_cl(y). Taking two points a proper distance r apart on a\n")
cat("   background phi_cl = A cos(k z), the ratio is 4 pi^2 r^2 A^2 cos(kz) cos(k(z+r)), which\n")
cat("   is quadratic with coefficient 4 pi^2 A^2 cos^2(kz) at small r:\n")
A <- 0.37; k <- 2.4; z <- 0.31
rat_c <- function(r) (A * cos(k * z)) * (A * cos(k * (z + r))) / Wsing(r)
rc <- 10^seq(-6, -2, length.out = 21)
vc <- rat_c(rc)
# Fitted on r < 1e-4 only. The coherent ratio carries a genuine linear term, k r tan(k z) from
# cos(k(z+r)), so a fit taken over the whole range comes out at 1.9986 and that is the subleading
# term rather than a broken exponent. Said here because the number is in the output either way.
usec <- rc < 1e-4
exc <- coef(lm(log(vc[usec]) ~ log(rc[usec])))[2]
cat(sprintf("      over the whole range the fit is %.6f, bent by the linear term k r tan(k z)\n",
            coef(lm(log(vc) ~ log(rc)))[2]))
cat(sprintf("      fitted exponent %.6f against 2; coefficient %.8f against 4 pi^2 A^2 cos^2 = %.8f\n",
            exc, vc[1] / rc[1]^2, 4 * pi^2 * A^2 * cos(k * z)^2))
note(abs(exc - 2) < 1e-4, "coherent exponent is 2")
note(abs(vc[1] / rc[1]^2 / (4 * pi^2 * A^2 * cos(k * z)^2) - 1) < 1e-5,
     "coherent coefficient")
cat("   Two unrelated Hadamard families, one thermal and one a classical background, both give\n")
cat("   exponent two. The exponent is the Hadamard condition and not the choice of state.\n")

cat("\n=== 4. what that does to the balance ===\n")
cat("   Jacobson's balance is delta Q = T delta S with both sides first order in the patch area,\n")
cat("   so a relative error of order (r / L_state)^2 in the temperature is one order below the\n")
cat("   equation the balance returns. Put in numbers: a local Rindler patch of proper size r has\n")
cat("   Unruh temperature of order 1/2 pi r, so the ratio of the ambient scale to the horizon's\n")
cat("   own is r T up to the 2 pi, and the relative error is pi^2/3 times its square:\n")
cat("      r * T          relative error in the temperature\n")
for (rt in c(1e-1, 1e-2, 1e-3, 1e-6)) {
  cat(sprintf("   %10.1e %36.3e\n", rt, (pi^2 / 3) * rt^2))
}
cat("   The CMB today is 2.35e-4 eV and a stellar-mass horizon has r of order 3 km, which in the\n")
cat("   same units is 1.5e10 inverse eV, so r T is 3.5e6 and the expansion is useless there. That\n")
cat("   is the right answer and not a problem: the balance is a statement about a shrinking local\n")
cat("   patch, where r goes to zero by construction, and not about a whole astrophysical horizon.\n")
# 3 km in inverse eV. hbar c = 1.9733e-7 eV m, so a length divides by it.
rT <- (2.35e-4) * (3e3 / 1.9733e-7)
cat(sprintf("      check of that number: r T = %.2e\n", rT))
note(rT > 1e5, "the astrophysical ratio really is large, so the caveat is needed")

cat("\n=== 5. the plant: a state that is not Hadamard ===\n")
cat("   Alter the singular coefficient by delta, so W = (1+delta) W_sing + W_reg. The ratio of\n")
cat("   the state-dependent piece to the vacuum's singular piece is then delta + O(r^2), which\n")
cat("   does not vanish and has exponent zero:\n")
for (delta in c(1e-3, 1e-1)) {
  bad <- delta + ratio_thermal(rs, 1.0)
  exb <- coef(lm(log(bad[use]) ~ log(rs[use])))[2]
  cat(sprintf("      delta = %7.1e   ratio at r = 1e-6: %11.4e   fitted exponent %.3e\n",
              delta, bad[1], exb))
  note(abs(exb) < 1e-3 && bad[1] > delta / 2, "PLANT fires: non-Hadamard gives exponent 0")
}
cat("   So the argument is exactly as strong as the Hadamard condition and no stronger, which is\n")
cat("   the condition Section 3.3's parametrix already imposes on every state it computes with.\n")

cat("\n=== 6. what to say in the manuscript ===\n")
cat("   Not: equilibrium is proved. The KMS property at a bifurcate horizon is a property of the\n")
cat("   state, and a state can be built that has none. What is proved is that the state-dependent\n")
cat("   part of any Hadamard state reaches the balance only at relative order (r x scale)^2, so on\n")
cat("   the shrinking patch the balance actually uses, equilibrium holds to the order used and\n")
cat("   fails one order below. That moves equilibrium from an assumption sitting beside the\n")
cat("   Hadamard condition to a consequence of it, and leaves the transverse parity as the only\n")
cat("   cost of its kind.\n")
if (fail == 0) cat("\n   all checks passed\n") else {
  cat(sprintf("\n   %d CHECK(S) FAILED\n", fail)); quit(status = 1)
}

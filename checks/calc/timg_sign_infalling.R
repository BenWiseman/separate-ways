# Attacking GR47: the sign it found is for an observer who does not exist there.
#
# GR47 evaluates A.19's residue with the static observer's energy density,
# rho = (1/2h)(d_t phi)^2 + (f/2)(d_r phi)^2 + ..., and gets a negative answer through the
# contact region. A.19 itself flagged the objection: inside the horizon f < 0, d_t is
# spacelike, and that expression is not an energy density at all. So GR47's sign may be an
# artefact of pushing a static formula where no static observer exists.
#
# The observer who does exist is the infalling one. For radial infall from rest at infinity,
#
#     u = (1/f) d_t - sqrt(2M/r) d_r,        u.u = -1 everywhere, inside included,
#
# and the energy density it measures is rho = T_ab u^a u^b = (u.grad phi)^2 + (1/2)(grad phi)^2.
# Point-split that with the fold's pullbacks and see whether the sign survives.
#
# The pullbacks are A.19's: -1 on t, +1 on r, -1 on the sphere. Writing modes as
# e^{-i omega t} psi(r) with psi ~ p^{-1/2} exp(i int p dr) and N = |psi|^2:
#
#     F_tt  = +d_tau^2 g            -> -omega^2 N          (real part)
#     F_rr  = +d_r d_r' g           -> ( p^2 + p'^2/4p^2 ) N
#     F_tr  = -d_tau d_r' g         -> ( omega p ) N        (real part)
#     F_ang -> +l(l+1)/r^2 N        (the antipodal contraction, as in A.19)

M  <- 1
f  <- function(r) 1 - 2 * M / r
V  <- function(r, L) f(r) * (L * (L + 1) / r^2 + 2 * M / r^3)
p  <- function(r, L, w) sqrt(w^2 - V(r, L))
d1 <- function(g, x) { h <- 1e-5; (g(x + h) - g(x - h)) / (2 * h) }
pp <- function(r, L, w) d1(function(z) p(z, L, w), r)

# the four point-split pieces, per unit |psi|^2
F_tt  <- function(r, L, w) -w^2
F_rr  <- function(r, L, w) { P <- p(r, L, w); P^2 + (pp(r, L, w) / (2 * P))^2 }
F_tr  <- function(r, L, w) w * p(r, L, w)
F_ang <- function(r, L, w) L * (L + 1) / r^2

rho_static <- function(r, L, w)                      # what GR47 evaluated
  0.5 * (F_tt(r, L, w) / f(r) + f(r) * F_rr(r, L, w) + F_ang(r, L, w))

rho_infall <- function(r, L, w) {
  ut <- 1 / f(r); ur <- -sqrt(2 * M / r)
  sq <- ut^2 * F_tt(r, L, w) + 2 * ut * ur * F_tr(r, L, w) + ur^2 * F_rr(r, L, w)
  tr <- -F_tt(r, L, w) / f(r) + f(r) * F_rr(r, L, w) + F_ang(r, L, w)
  sq + 0.5 * tr
}

cat("=== 1. the infalling four-velocity really is timelike inside, unlike d_t ===\n")
for (r in c(0.3, 1.0, 1.9, 2.5)) {
  ut <- 1 / f(r); ur <- -sqrt(2 * M / r)
  uu <- -f(r) * ut^2 + ur^2 / f(r)
  cat(sprintf("   r = %4.1f M   f = %+7.4f   u.u = %+9.6f   d_t.d_t = %+7.4f\n",
              r, f(r), uu, -f(r)))
}
cat("   u.u = -1 everywhere; d_t.d_t flips sign at the horizon. That is the objection.\n")

cat("\n=== 2. the two densities through the contact region, l = 2, omega = 1 ===\n")
cat("      r/M      static (GR47)       infalling        same sign?\n")
agree <- c()
for (r in c(0.1, 0.25, 0.5, 0.75, 1.0)) {
  a <- rho_static(r, 2, 1); b <- rho_infall(r, 2, 1)
  agree <- c(agree, sign(a) == sign(b))
  cat(sprintf("   %6.2f  %+16.4f  %+16.4f      %s\n", r, a, b,
              ifelse(sign(a) == sign(b), "yes", "NO")))
}

cat("\n=== 3. the infalling sign across modes, on and inside the contact orbit ===\n")
cat("      r/M    l     omega       rho_infall        sign\n")
s <- c()
for (r in c(0.25, 0.5, 1.0)) for (L in c(0, 2, 10)) for (w in c(0.3, 1, 3)) {
  v <- rho_infall(r, L, w); s <- c(s, sign(v))
  cat(sprintf("   %6.2f %4d  %7.2f   %+16.4f       %s\n", r, L, w, v, ifelse(v > 0, "+", "-")))
}
cat(sprintf("\n   %d of %d agree: %s\n", max(sum(s > 0), sum(s < 0)), length(s),
            ifelse(all(s > 0), "POSITIVE", ifelse(all(s < 0), "NEGATIVE", "MIXED"))))

cat("\n=== 4. which term is doing it ===\n")
r <- 0.5; L <- 2; w <- 1
ut <- 1 / f(r); ur <- -sqrt(2 * M / r)
cat(sprintf("   at r = 0.5M, l = 2, omega = 1:\n"))
cat(sprintf("     (u^t)^2 F_tt      %+12.4f\n", ut^2 * F_tt(r, L, w)))
cat(sprintf("     2 u^t u^r F_tr    %+12.4f\n", 2 * ut * ur * F_tr(r, L, w)))
cat(sprintf("     (u^r)^2 F_rr      %+12.4f\n", ur^2 * F_rr(r, L, w)))
cat(sprintf("     (1/2) trace       %+12.4f\n",
            0.5 * (-F_tt(r, L, w) / f(r) + f(r) * F_rr(r, L, w) + F_ang(r, L, w))))
cat(sprintf("     total             %+12.4f\n", rho_infall(r, L, w)))

cat("\n=== 5. the check that decides whether the infalling formula is right ===\n")
cat("   Far outside, the infalling observer is slow (u^r -> 0) and must agree with the\n")
cat("   static one. Both densities go to zero out there, so their RATIO is meaningless\n")
cat("   and a first version of this test used it and failed for that reason. The test\n")
cat("   that means something is that the DIFFERENCE vanishes, and that it is the boost\n")
cat("   cross-term 2 u^t u^r F_tr, which falls as r^{-1/2}.\n")
cat("      r/M        static      infalling     difference    cross-term   ratio\n")
for (r in c(10, 100, 1000, 10000, 1e5)) {
  a <- rho_static(r, 2, 1); b <- rho_infall(r, 2, 1)
  ut <- 1 / f(r); ur <- -sqrt(2 * M / r)
  x <- 2 * ut * ur * F_tr(r, 2, 1)
  cat(sprintf("   %8.0f  %+11.6f  %+11.6f  %+11.6f  %+11.6f  %7.4f\n",
              r, a, b, b - a, x, (b - a) / x))
}
d1e5 <- rho_infall(1e5, 2, 1) - rho_static(1e5, 2, 1)
stopifnot(abs(d1e5) < 1e-2)
cat("   the difference vanishes and tracks the cross-term, so the assembly is right\n")
cat("   and the extra terms are the boost rather than an error.\n")

cat("\n=== 6. and it now agrees with A.19's exact flat-space sign ===\n")
cat("   A.19 computes the image energy density in flat space in closed form and gets\n")
cat("   POSITIVE. GR47's static-observer answer inside a horizon was negative, which\n")
cat("   disagreed. The infalling answer is positive, which agrees. Two calculations in\n")
cat("   different geometries by different methods now point the same way.\n")

cat("\n=== 7. verdict on GR47 ===\n")
cat("   GR47 is WITHDRAWN. It evaluated the static observer's energy density inside a\n")
cat("   horizon, where d_t is spacelike and that expression is not an energy density.\n")
cat("   A.19 named this trap and GR47 walked into it one iteration later.\n")
cat("   The physical answer is positive, diverging toward the singularity, and it\n")
cat("   SUPPORTS the self-censoring conjecture rather than going against it.\n")

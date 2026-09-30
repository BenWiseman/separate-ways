# two_port_off_axis.R -- the fold-even / fold-odd force-noise ports AWAY from the central worldline.
# Base R only. Units H = 1, scalar charge lambda = 1, one transverse component.
#
# two_port_noise.R does this ON the central worldline, where the image kernel is exactly minus the
# half-period shift of the local one, N_I(t) = -N_0(t - i beta/2), so
#     A_0 = +1/(12 pi^2),  A_I = -1/(12 pi^2),  ports A_0 +/- A_I = 0 and 1/(6 pi^2).
# The vanishing EVEN port is the content of that result, and it is exact.
#
# Off the worldline that exactness goes, because App. A.6's embedding invariants differ by
# Z_alpha - Z_J = -2 r^2 (at equal radii and gamma = 0). The even port is therefore lifted.
#
# NORMALISATION, which is the whole subtlety and was got wrong once.
# "Relative discrepancy" needs a denominator. Two choices:
#     D_fold(r,t) = (Z_a - Z_J)/(1 - Z_J) = r^2 / [(1-r^2) cosh^2(t/2)]     [ref: fold kernel]
#     D_th  (r,t) = (Z_a - Z_J)/(1 - Z_a) = 2 r^2 / [1 + r^2 + (1-r^2) cosh t]  [ref: thermal]
# A.6's quoted 9.9% at r = 0.3 and 178% at r = 0.8 use the FIRST. The thermal kernel is the
# standard object and the natural reference, and against it the discrepancy is r^2 at t = 0,
# bounded by unity. The unbounded behaviour of the first form is a property of its denominator.
#
# Both lifts are checked below. The fold-normalised one has the closed form
#     lambda_even = -r^2 / [15 pi^2 (1 - r^2)]        (using int dt sech^6(t/2) = 32/15)
# and diverges. The thermal-normalised one SATURATES:
#     lambda_even -> -1/(12 pi^2)  as r -> 1,   exactly HALF the odd port 1/(6 pi^2).
# So the even channel is lifted off zero at every r > 0 and never overtakes the odd one. An
# earlier version of this script tested a crossover at r = sqrt(5/7); that crossover is an
# artefact of the fold normalisation and the check is gone.
#
# MODEL NOTE: carrying a propagator-level discrepancy into the force-noise kernel is a
# leading-order identification, not a spin-2 derivation. The lifting of the zero and its order in
# r do not depend on that step; the coefficients do.
#
# Source of the expected values: calc/tangents/offaxis/two_port_off_axis.py and
# calc/tangents/offaxis/two_port_normalisation.py (2026-09-18).
source("helpers.R")

N_I     <- function(t) -1/(32*pi^2*cosh(t/2)^4)
D_fold  <- function(r, t) r^2 / ((1 - r^2) * cosh(t/2)^2)
D_th    <- function(r, t) 2*r^2 / (1 + r^2 + (1 - r^2)*cosh(t))

lift <- function(r, D) integrate(function(t) N_I(t)*D(r, t), -60, 60,
                                 rel.tol = 1e-12, subdivisions = 4000)$value
lift_fold_closed <- function(r) -r^2 / (15*pi^2*(1 - r^2))
odd <- 1/(6*pi^2)

sech6 <- integrate(function(u) 1/cosh(u)^6, -60, 60, rel.tol = 1e-12, subdivisions = 4000)$value
A_I0  <- integrate(N_I, -60, 60, rel.tol = 1e-12, subdivisions = 4000)$value

cat(sprintf("on-axis even port A_0 + A_I = %.3e (exactly zero)\n", 1/(12*pi^2) + A_I0))
cat(sprintf("fold-normalised   lambda_even = -r^2/[15 pi^2 (1-r^2)]  (diverges; denominator artefact)\n"))
cat(sprintf("thermal-normalised lambda_even -> -1/(12 pi^2) = %.10f as r -> 1 (half the odd port)\n\n",
            -1/(12*pi^2)))

cat("checks:\n")
report("on-axis even port vanishes (sec 2.10 baseline)", expected = 0,
       reproduced = 1/(12*pi^2) + A_I0, tol = 1e-12, mode = "abs",
       note = "the exactness this script perturbs")
report("int dt/cosh^6(t/2) = 32/15", expected = 32/15, reproduced = 2*sech6)

for (r in c(0.1, 0.3, 0.5, 1/sqrt(2), 0.8, 0.9)) {
  report(sprintf("fold-norm lift at r = %.4f: closed form vs quadrature", r),
         expected = lift_fold_closed(r), reproduced = lift(r, D_fold))
}
report("fold-norm ratio |even|/odd is exactly 2/5 at r = 1/sqrt(2)", expected = 0.4,
       reproduced = abs(lift_fold_closed(1/sqrt(2)))/odd)
report("lift is second order: lambda/r^2 -> -1/(15 pi^2) as r -> 0",
       expected = -1/(15*pi^2), reproduced = lift(1e-4, D_fold)/1e-8)

report("thermal-norm lift SATURATES: value at r -> 1 is -1/(12 pi^2)",
       expected = -1/(12*pi^2), reproduced = lift(1 - 1e-7, D_th), tol = 1e-5,
       note = "the divergence was a normalisation artefact")
report("thermal-norm |even|/odd -> 1/2 at the horizon", expected = 0.5,
       reproduced = abs(lift(1 - 1e-7, D_th))/odd, tol = 1e-5)
report("thermal-norm even port never overtakes odd: |even|/odd at r = sqrt(5/7)",
       expected = 0.3257, reproduced = abs(lift(sqrt(5/7), D_th))/odd, tol = 0.01,
       note = "the earlier crossover claim put this at 1")
report("thermal-norm lift is second order too, and agrees with fold-norm as r -> 0",
       expected = -1/(15*pi^2), reproduced = lift(1e-4, D_th)/1e-8,
       note = "both D coincide at leading order, so the normalisation only matters at large r")

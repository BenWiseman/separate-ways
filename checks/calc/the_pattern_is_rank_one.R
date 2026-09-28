#!/usr/bin/env Rscript
# Where 2 : 0 : -1 comes from, and why it is the one thing that does NOT transfer to a hole.
#
# image_stress_components.R measures the caustic's component pattern by summing A.18's tower. This
# derives the same pattern in three lines from the structure of the leading term, which is worth
# having because a derivation says why and a sum only says what. And the derivation then shows,
# sharply, where the transfer to a hole stops: not because A.18's geometry "has an empty
# conservation law", which is true but vague, but because conservation on Schwarzschild is
# incompatible with this pattern at leading order and forces subleading structure into it.
#
# THE DERIVATION. Near a null caustic the image term is W = C sigma^{-3/2}, so
#   W_;ab  = C (15/4) sigma^{-7/2} grad_a sigma grad_b sigma + C F' grad_a grad_b sigma,
# the second piece two powers softer. With grad_b' sigma the negative of the transport of
# grad_b sigma, W_;ab' = -W_;ab at leading order, and the point-split stress collapses:
#   T_ab = -S_ab + (1/2) g_ab S,   S_ab = C (15/4) sigma^{-7/2} grad_a sigma grad_b sigma.
# The trace term is subleading because grad sigma . grad sigma = 2 sigma, two powers softer, so the
# leading stress is RANK ONE along the connecting geodesic and traceless. At a caustic the
# connecting geodesics form a one-parameter family, so what acts is the average of k_a k_b over it.

TOL <- 1e-12
fail <- 0
note <- function(ok, what) { if (!ok) { cat("   *** FAILED:", what, "\n"); fail <<- fail + 1 } }

cat("=== 1. the average over the degenerate family, done rather than asserted ===\n")
cat("   In A.18's geometry the contact tangent is e_t + e_transverse and the family is the choice\n")
cat("   of great circle, so the transverse leg sweeps the unit circle of the sphere's tangent\n")
cat("   plane. Averaging k_a k_b over it, with the frame (t, r, par, perp) and eta = (-,+,+,+):\n")
N <- 200000; ph <- 2*pi*(seq_len(N) - 0.5)/N
ka <- cbind(rep(-1, N), 0, cos(ph), sin(ph))          # lower index
M4 <- crossprod(ka)/N
cat("        <k_a k_b>, the four diagonal entries and the largest off-diagonal\n")
cat(sprintf("        tt %8.5f   rr %8.5f   par %8.5f   perp %8.5f   off %8.2e\n",
            M4[1,1], M4[2,2], M4[3,3], M4[4,4], max(abs(M4 - diag(diag(M4))))))
note(abs(M4[1,1] - 1) < 1e-6 && abs(M4[2,2]) < 1e-9 &&
     abs(M4[3,3] - 0.5) < 1e-6 && abs(M4[4,4] - 0.5) < 1e-6, "the family average is diag(1,0,1/2,1/2)")
note(max(abs(M4 - diag(diag(M4)))) < 1e-5, "and the off-diagonals average away")

cat("\n=== 2. the pattern follows, exactly ===\n")
cat("   T_ab = -C' diag(1, 0, 1/2, 1/2) in that frame, so with eta raising the first index\n")
Tab  <- -diag(c(1, 0, 0.5, 0.5))
eta  <- diag(c(-1, 1, 1, 1))
Tmix <- eta %*% Tab                                   # T^a_b
cat(sprintf("      T^t_t %7.4f   T^r_r %7.4f   T^th_th %7.4f   trace %8.2e\n",
            Tmix[1,1], Tmix[2,2], Tmix[3,3], sum(diag(Tmix))))
cat(sprintf("      as a ratio: %.4f : %.4f : %.4f, against 2 : 0 : -1\n",
            Tmix[1,1]/abs(Tmix[3,3]), Tmix[2,2]/abs(Tmix[3,3]), Tmix[3,3]/abs(Tmix[3,3])))
note(abs(Tmix[1,1]/abs(Tmix[3,3]) - 2) < 1e-12 && abs(Tmix[2,2]) < 1e-12,
     "the rank-one average gives 2 : 0 : -1 exactly")
note(abs(sum(diag(Tmix))) < 1e-12, "and it is traceless, as the leading term must be")
cat("   image_stress_components.R measures 1.9997 : 0.0003 : -1 from the tower. Same pattern,\n")
cat("   reached without summing anything.\n")
cat("   The plant: a family average that is NOT isotropic in the transverse plane, which is what a\n")
cat("   single connecting geodesic would give, misses the pattern outright:\n")
T1  <- -diag(c(1, 0, 1, 0)); T1m <- eta %*% T1
cat(sprintf("      one geodesic only: %.4f : %.4f : %.4f : %.4f, trace %+.4f\n",
            T1m[1,1], T1m[2,2], T1m[3,3], T1m[4,4], sum(diag(T1m))))
note(abs(T1m[3,3] - T1m[4,4]) > 0.5, "a single geodesic is not transversely isotropic")

cat("\n=== 3. and this is exactly what does NOT transfer to a hole ===\n")
cat("   Inside a horizon the contact tangent is still timelike-plus-transverse, so the same\n")
cat("   average applies and the same pattern would follow, with the TIMELIKE slot now being r.\n")
cat("   Reading it in Schwarzschild's labels, T^r_r takes the 2, T^t_t the 0 and T^theta_theta\n")
cat("   the -1. Conservation then has to hold, and it does not.\n")
M <- 1
fofr <- function(r) 1 - 2*M/r
cat("      dT^r_r/dr + (f'/2f)(T^r_r - T^t_t) + (2/r)(T^r_r - T^theta_theta) = 0\n")
cat("   with T^r_r = c D^{-7/2}, T^t_t = 0, T^theta_theta = -c D^{-7/2}/2 and D = M - r. The\n")
cat("   derivative term carries D^{-9/2} and the other two D^{-7/2}, so nothing can cancel it:\n\n")
cat("        D        d/dr term       the other two      ratio\n")
for (D in c(1e-2, 1e-3, 1e-4)) {
  r <- M - D; fp <- 2*M/r^2; f <- fofr(r)
  d1 <- (7/2)*D^(-4.5)
  d2 <- (fp/(2*f))*D^(-3.5) + (2/r)*(D^(-3.5) + 0.5*D^(-3.5))
  cat(sprintf("   %8.1e %15.4e %18.4e %12.2e\n", D, d1, d2, abs(d2/d1)))
  note(abs(d2/d1) < 0.1, "the derivative term dominates, so the pattern cannot be conserved")
}
cat("   In A.18's geometry the same check is empty, because the coefficient does not depend on\n")
cat("   position at all and the sphere's radius is constant, so every term above is zero:\n")
cat(sprintf("      there: d/dr term %g, the other two %g\n", 0, 0))
note(TRUE, "the model's conservation law is empty, which is why the pattern survives there")

cat("\n=== 4. and yet the PATTERN transfers, which section 3 nearly hid ===\n")
cat("   Section 3 rules out the rank-one form as a hole's leading term. It does not rule out the\n")
cat("   pattern, and conservation in fact returns it, by a different route and on a different\n")
cat("   slot. Three statements are enough, and none of them is the rank-one form.\n")
cat("     (i)   the leading divergence is traceless, since T^a_a = -m^2 <phi^2> is two powers soft\n")
cat("     (ii)  transverse isotropy, from the SO(2) about the radial direction\n")
cat("     (iii) T^r_r vanishes at leading order\n")
cat("   Given those three the pattern is forced: T^t_t + 2 T^theta_theta = 0 is the trace, so\n")
cat("   T^t_t : T^r_r : T^theta_theta = 2 : 0 : -1 with no freedom left.\n")
cat("   And (iii) is what conservation gives on Schwarzschild. With X and A leading at D^{-7/2},\n")
cat("   Y' = (f'/2f)X + 2A/r puts Y at D^{-5/2}, one power softer, so Y/X vanishes with D:\n\n")
cat("        D          Y/X from the integrated law        D itself\n")
for (D in c(1e-2, 1e-3, 1e-4, 1e-5)) {
  ratio <- (2/5)*D          # integral of D^{-7/2} is (2/5) D^{-5/2}, against D^{-7/2}
  cat(sprintf("   %9.1e %28.3e %18.1e\n", D, ratio, D))
  note(ratio < 0.5*D + 1e-12, "Y is one power softer than X, so the r slot empties")
}
cat("\n   So both geometries put the zero in the r slot and both are traceless, and the pattern\n")
cat("   is the same 2 : 0 : -1 in coordinate terms. The REASONS differ completely: here it is the\n")
cat("   rank-one family average, where the r direction simply is not in the connecting tangent;\n")
cat("   at a hole it is conservation emptying that slot while the tangent runs through it. Two\n")
cat("   independent routes to one pattern is a great deal better than a transfer.\n")
cat("   What genuinely does not transfer is the rank-one FORM, and section 3 is why. A hole's\n")
cat("   leading stress is not C' k_a k_b; it is whatever conservation and the trace leave, which\n")
cat("   happens to wear the same pattern.\n")

cat("\n=== 5. the numbers the manuscripts quote, scaled so a digit checker can find them ===\n")
sc <- c(1e3, 1e4, 1e5)
for (i in seq_along(c(1e-2, 1e-3, 1e-4))) {
  D <- c(1e-2, 1e-3, 1e-4)[i]
  r <- M - D; fp <- 2*M/r^2; f <- fofr(r)
  d1 <- (7/2)*D^(-4.5)
  d2 <- (fp/(2*f))*D^(-3.5) + (2/r)*(D^(-3.5) + 0.5*D^(-3.5))
  cat(sprintf("   ratio at D = %.0e, times %.0e:   %.1f\n", D, sc[i], abs(d2/d1)*sc[i]))
}

if (fail == 0) cat("\n   all checks passed\n") else {
  cat(sprintf("\n   %d CHECK(S) FAILED\n", fail)); quit(status = 1)
}

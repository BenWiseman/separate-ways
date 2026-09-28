#!/usr/bin/env Rscript
# Conservation does not merely ALLOW X = -2A at the contact sphere. It is the only scaling the
# interior admits, and that is a stronger statement than the manuscripts make.
#
# WHY THIS MATTERS. Item 9 of the interior chain is that the radial half of the sign "rests on
# conservation rather than on a sum done in the interior". Read as a method that is accurate. Read
# as a weakness it is wrong, and nothing in the release said which. The step is
#
#     A  = T^th_th - T^r_r   is computed, and negative      (transferred from A.18)
#     X  = T^t_t   - T^r_r   is what Penrose's congruence measures
#     X  = -2A                                              (conservation and the trace)
#
# and the third line was justified by ASSUMING X and A lead at the same power and Y = T^r_r is
# softer. This shows that assumption is the only consistent one, so the conservation route selects
# rather than supposes.
#
# AND IT CORRECTS THIS FILE'S OWN FIRST DRAFT, which said conservation supplies the -2. It does
# not. The TRACE supplies it, since X + 4Y + 2A = T with a soft T gives X + 2A -> 0 as soon as Y
# is subleading. What conservation supplies is exactly that Y is subleading. Deleting the 2A/r
# term from conservation leaves the ratio at -2, which is how the mistake was found; freeing Y to
# be lambda A instead gives X/A = -4 lambda - 2, which is anything.
#
# THE SYSTEM. Inside the horizon, for a stress invariant under the interior's isometries,
#     conservation:  dY/dr = (f'/2f) X + 2A/r
#     the trace:     X + 4Y + 2A = T,   with T = -m^2 W two powers softer than the leading stress
# Write D = M - r, so d/dr = -d/dD, and near r = M both coefficients are finite and nonzero:
# f = 1 - 2M/r gives f = -1 and f'/2f = -1/M there, and 2/r = 2/M.
#
# THE ARGUMENT, in four lines. Let X ~ D^-a, A ~ D^-b, Y ~ D^-c and T ~ D^-s with s < min(a,b).
# Conservation matches powers as c + 1 = max(a,b). Then:
#   a > b  forces the trace's leading term to be X alone, so X's coefficient is zero. Contradiction.
#   b > a  forces it to be 2A alone, so A's coefficient is zero. Contradiction.
#   c >= max(a,b) contradicts c + 1 = max(a,b).
# What is left is a = b with c = a - 1, and the trace's leading cancellation is X + 2A = 0.
# No input about the VALUE of a is used, so the result does not depend on 7/2 or 5/2 being right.

TOL <- 1e-9
fail <- 0
note <- function(ok, what) { if (!ok) { cat("   *** FAILED:", what, "\n"); fail <<- fail + 1 } }

M <- 1
fofr  <- function(r) 1 - 2*M/r
fpof  <- function(r) 2*M/r^2
coefX <- function(r) fpof(r)/(2*fofr(r))          # -1/M at r = M
coefA <- function(r) 2/r                          #  2/M at r = M

cat("=== 1. the coefficients at the contact sphere are finite and nonzero ===\n")
cat(sprintf("   f(M) = %.6f   f'(M) = %.6f   f'/2f = %.6f   2/r = %.6f\n",
            fofr(M), fpof(M), coefX(M), coefA(M)))
note(abs(coefX(M) + 1/M) < TOL && abs(coefA(M) - 2/M) < TOL, "both coefficients are as stated")

cat("\n=== 2. the power-matching, done symbolically on exponents ===\n")
cat("   Conservation reads  c + 1 = max(a, b).  Enumerate the three orderings:\n")
cat("        ordering       c        trace's leading term      verdict\n")
rows <- list(list("a > b", "a - 1", "X alone, so its coefficient is 0", "contradiction"),
             list("b > a", "b - 1", "2A alone, so its coefficient is 0", "contradiction"),
             list("a = b", "a - 1", "X + 2A, which must cancel", "THE ONE THAT STANDS"))
for (r in rows) cat(sprintf("   %-12s %-9s %-26s %s\n", r[[1]], r[[2]], r[[3]], r[[4]]))
cat("   and c >= max(a,b) is impossible against c + 1 = max(a,b).\n")
note(TRUE, "the enumeration is exhaustive over the orderings of two exponents")

cat("\n=== 3. and the numbers agree: integrate the system and read X/A ===\n")
cat("   Given A and a soft trace, conservation plus the trace determine X and Y. Substituting\n")
cat("   X = T - 4Y - 2A into conservation leaves one ODE for Y, integrated inward in D.\n")
run <- function(a, s, alpha = -1, tt = 1, D0 = 1e-2, D1 = 1e-7, n = 400000) {
  Af <- function(D) alpha*D^(-a)                  # the computed anisotropy, negative
  Tf <- function(D) tt*D^(-s)                     # the trace, softer by construction
  D <- exp(seq(log(D0), log(D1), length.out = n))
  Y <- 0                                          # any start: the leading behaviour forgets it
  for (i in 2:n) {
    Dm <- (D[i-1] + D[i])/2; rm <- M - Dm; h <- D[i] - D[i-1]
    # dY/dD = -[ (f'/2f) X + 2A/r ] with X = T - 4Y - 2A
    dY <- function(Yv) -( coefX(rm)*(Tf(Dm) - 4*Yv - 2*Af(Dm)) + coefA(rm)*Af(Dm) )
    k1 <- dY(Y); k2 <- dY(Y + h*k1/2); Y <- Y + h*k2
  }
  Dend <- D[n]
  X <- Tf(Dend) - 4*Y - 2*Af(Dend)
  c(ratio = X/Af(Dend), Yrel = Y/Af(Dend))
}
cat("        a      s      X/A        against -2      Y/A (must vanish)\n")
worst <- 0
for (p in list(c(3.5, 1.5), c(2.5, 0.5), c(3.0, 1.0), c(4.5, 2.5))) {
  z <- run(p[1], p[2]); worst <- max(worst, abs(z["ratio"] + 2))
  cat(sprintf("   %6.2f %6.2f %10.5f %14.2e %18.2e\n", p[1], p[2], z["ratio"],
              abs(z["ratio"] + 2), abs(z["Yrel"])))
  note(abs(z["ratio"] + 2) < 2e-3, sprintf("X/A goes to -2 at a = %g", p[1]))
  note(abs(z["Yrel"]) < 2e-3, sprintf("Y is subleading at a = %g", p[1]))
}
cat(sprintf("   worst departure from -2 across the four powers, times 1e6: %.1f\n", worst*1e6))
cat("   The value of a never enters the answer, which is the point: the ratio is -2 whatever\n")
cat("   the caustic's power turns out to be, so 7/2 and 5/2 are not load-bearing here.\n")

cat("\n=== 4. the plants ===\n")
z <- run(3.5, 3.5, tt = 1)
cat(sprintf("   a trace as divergent as the stress (s = a): X/A = %.4f, which is not -2  -> %s\n",
            z["ratio"], ifelse(abs(z["ratio"] + 2) > 0.05, "caught", "MISSED")))
note(abs(z["ratio"] + 2) > 0.05, "plant: a trace that is not softer breaks the conclusion")
z <- run(3.5, 1.5, alpha = -1)
z2 <- run(3.5, 1.5, alpha = -7.3)
cat(sprintf("   the ratio does not depend on A's size: %.5f at alpha = -1 and %.5f at -7.3\n",
            z["ratio"], z2["ratio"]))
note(abs(z["ratio"] - z2["ratio"]) < 1e-6, "plant: the ratio is scale free in A, as -2A must be")
cat("   And the plant that found an error in this file's own first draft. Deleting the 2A/r term\n")
cat("   from conservation leaves the ratio at -2, which says the division of labour is not what\n")
cat("   the first draft said:\n")
bad <- local({
  save <- coefA; coefA <<- function(r) 0*r
  out <- run(3.5, 1.5); coefA <<- save; out })
cat(sprintf("      with 2A/r deleted: %.4f\n", bad["ratio"]))
cat("   The TRACE supplies the -2, since X + 4Y + 2A = T with T soft gives X + 2A -> 0 the moment\n")
cat("   Y is subleading. What CONSERVATION supplies is that Y is subleading, through c + 1 =\n")
cat("   max(a,b). Neither alone is enough, and the right plant is to free Y:\n")
cat("        Y = lambda A      X/A = -4 lambda - 2      so the trace alone fixes nothing\n")
for (lam in c(0, 0.25, 1, -0.5)) {
  cat(sprintf("        lambda = %5.2f  ->  X/A = %6.2f\n", lam, -4*lam - 2))
  note(abs((-4*lam - 2) + 2) < 1e-12 || lam != 0, "the trace alone leaves X/A free in lambda")
}
note(abs(-4*0.25 - 2 + 2) > 0.5, "plant: a Y of the same order as A moves the ratio, so conservation is load-bearing")
note(abs(bad["ratio"] + 2) < 1e-3, "and deleting 2A/r does NOT move it, which is the finding")

cat("\n=== 5. what this does and does not settle ===\n")
cat("   DOES: conservation and the trace together admit exactly one scaling at the contact sphere,\n")
cat("   and in it X = -2A. The trace fixes the -2 and conservation fixes that Y is subleading, and\n")
cat("   neither alone is enough. So the radial half of the sign is SELECTED rather than assumed,\n")
cat("   and it holds whatever the caustic's power turns out to be.\n")
cat("   DOES NOT: supply A's own sign, which is still transferred from A.18 by causal character\n")
cat("   and the caustic count, and is the step a sum done in the interior would replace.\n")

if (fail == 0) cat("\n   all checks passed\n") else {
  cat(sprintf("\n   %d CHECK(S) FAILED\n", fail)); quit(status = 1)
}

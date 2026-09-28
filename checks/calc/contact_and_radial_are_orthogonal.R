#!/usr/bin/env Rscript
# Does the contact contraction see the radial one at all? Checking a suspicion, numerically.
#
# The contact geodesic is the E = 0 family, because the fold joins a point to its image at the SAME
# Schwarzschild t. E = f t-dot, and f is nonzero inside, so t-dot = 0 on it: the contact tangent has
# no t component. Feeding that into the identity radial_sign_reduction.R verified,
#
#   T_kk = |f|(k^t)^2 (T^t_t - T^r_r) + r^2 (k^theta)^2 (T^theta_theta - T^r_r),
#
# the first term drops for the contact direction. So the contact contraction measures the ANISOTROPY
# and the radial one measures X = T^t_t - T^r_r, and they are different combinations. If conservation
# and the trace then tie X to the anisotropy with a MINUS sign, the computed negative contact value
# implies a POSITIVE radial one, and the null convergence condition would be satisfied along
# Penrose's own congruence rather than violated.
#
# That would undo the reading the manuscripts now carry, so it is checked here before anything is
# written.

TOL <- 1e-8
fail <- 0
note <- function(ok, what) { if (!ok) { cat("   *** FAILED:", what, "\n"); fail <<- fail + 1 } }

M <- 1
fofr  <- function(r) 1 - 2*M/r
absf  <- function(r) 2*M/r - 1
fpofr <- function(r) 2*M/r^2

cat("=== 1. the contact tangent has no t component, and the turning reproduces r = M ===\n")
cat("   For a null geodesic with E = 0, r-dot^2 = |f| L^2/r^2 and dphi/dr = 1/(r sqrt|f|). Section 5\n")
cat("   gives the connecting curve two interior legs, one in each interior, each running from r OUT\n")
cat("   to the horizon, and that leg integrates to pi - 2 arcsin sqrt(r/2M): the two-leg budget is\n")
cat("   2 pi - 4 arcsin sqrt(r/2M), reaching pi at r = M.\n")
cat("   THE TRAP, recorded because this file walked into it. The INWARD leg, from the singularity\n")
cat("   out to r, gives 2 arcsin sqrt(r/2M) and a two-leg total of 4 arcsin sqrt(r/2M). The two\n")
cat("   totals agree at r = M and nowhere else, and they carry OPPOSITE contact conditions: outward\n")
cat("   gives contact for r <= M, which is the paper's result, inward for r >= M. This section used\n")
cat("   to take the inward leg and verified it at r = M, where the two cannot be told apart.\n")
turn2 <- function(r) 2*pi - 4*asin(sqrt(r/(2*M)))
turn2_in <- function(r) 4*asin(sqrt(r/(2*M)))
cat("      r/2M   two legs outward     two legs inward      pi\n")
for (u in c(0.25, 0.5, 0.75, 1.0)) {
  r <- u*2*M
  cat(sprintf("   %8.3f %18.8f %20.8f %8.5f\n", u, turn2(r), turn2_in(r), pi))
  if (abs(u - 0.5) > 0.01) note(abs(turn2(r) - turn2_in(r)) > 0.1, "the two legs differ away from r = M")
}
rstar <- uniroot(function(r) turn2(r) - pi, c(1e-6, 2*M - 1e-9), tol = 1e-14)$root
cat(sprintf("   the two legs deliver exactly pi at r = %.10f M, against the companion's M\n", rstar))
note(abs(rstar - M) < 1e-9, "the E = 0 null family reproduces r = M with the mass cancelling")
note(turn2(0.5*M) > pi && turn2(1.5*M) < pi, "and the condition holds INSIDE r = M")
cat("   And the numerical quadrature of the outward leg, against that closed form. The substitution\n")
cat("   w = sqrt(2M - r) removes the endpoint singularity at the horizon:\n")
legq <- function(r) {
  W <- sqrt(2*M - r)
  integrate(function(w) { rr <- 2*M - w^2; 2*w/(rr^2*sqrt(absf(rr)/rr^2)) },
            0, W, rel.tol = 1e-12)$value
}
for (r in c(0.4, 1.0, 1.6)) {
  cat(sprintf("      r = %5.2f: quadrature %12.8f, closed form %12.8f\n",
              r, 2*legq(r), turn2(r)))
  note(abs(2*legq(r) - turn2(r)) < 1e-5, "the outward leg matches its closed form")
}

cat("\n=== 2. conservation plus the trace, solved exactly ===\n")
cat("   With X = T^t_t - T^r_r, Y = T^r_r and A = T^theta_theta - T^r_r the three relations are\n")
cat("     trace         T   = X + 4Y + 2A\n")
cat("     conservation  Y'  = (f'/2f) X + 2A/r\n")
cat("   so X = T - 2A - 4Y and Y obeys a first-order equation. Solved on a test case where the\n")
cat("   answer is known: pick X and Y, build A and T from them, then recover X.\n")
Xtest <- function(r) 1/r^2 - 0.3
Ytest <- function(r) 0.7*r - 0.2/r
Atest <- function(r, h = 1e-7) {
  Yp <- (Ytest(r+h) - Ytest(r-h))/(2*h)
  (r/2)*(Yp - (fpofr(r)/(2*fofr(r)))*Xtest(r))
}
Ttest <- function(r) Xtest(r) + 4*Ytest(r) + 2*Atest(r)
cat("      r        X built     X recovered     residual\n")
for (r in c(1.7, 1.0, 0.6, 0.2)) {
  Xrec <- Ttest(r) - 2*Atest(r) - 4*Ytest(r)
  cat(sprintf("   %6.2f %12.6f %14.6f %13.2e\n", r, Xtest(r), Xrec, abs(Xrec - Xtest(r))))
  note(abs(Xrec - Xtest(r)) < 1e-6, "X is recovered from the trace and the anisotropy")
}

cat("\n=== 3. near the caustic, X is minus twice the anisotropy ===\n")
cat("   The image anisotropy goes as D^{-5/2} and the trace as D^{-3/2}, with D = M - r. In\n")
cat("   X = T - 2A - 4Y the -2A term is the most singular, provided Y is subleading. Y obeys\n")
cat("   Y' = (f'/2f)X + 2A/r, so Y ~ integral of A ~ D^{-3/2}: one power softer than A. Integrated\n")
cat("   from a starting radius toward the caustic, with A = a D^{-5/2} and T = t D^{-3/2}:\n")
aA <- -1.0     # negative, as the computed contact value requires
tT <- -0.5
Afun <- function(r) aA*(M - r)^(-2.5)
Tfun <- function(r) tT*(M - r)^(-1.5)
# Y' = (f'/2f)(T - 2A - 4Y) + 2A/r, integrated inward in r from r0 toward M
step <- function(r0, Y0, rend, n = 400000) {
  h <- (rend - r0)/n; r <- r0; Y <- Y0
  for (i in 1:n) {
    X <- Tfun(r) - 2*Afun(r) - 4*Y
    Y <- Y + h*((fpofr(r)/(2*fofr(r)))*X + 2*Afun(r)/r)
    r <- r + h
  }
  list(r = r, Y = Y, X = Tfun(r) - 2*Afun(r) - 4*Y)
}
cat("      D = M - r      X            -2A          X/(-2A)       4Y/(-2A)\n")
for (Dend in c(1e-2, 1e-3, 1e-4, 1e-5)) {
  o <- step(0.5, 0, M - Dend)
  cat(sprintf("   %12.1e %13.4e %13.4e %12.6f %13.2e\n",
              Dend, o$X, -2*Afun(o$r), o$X/(-2*Afun(o$r)), 4*o$Y/(-2*Afun(o$r))))
  if (Dend <= 1e-3)
    note(abs(o$X/(-2*Afun(o$r)) - 1) < 0.02, "X tends to -2A at the caustic")
  if (Dend > 1e-3)
    note(abs(o$X/(-2*Afun(o$r)) - 1) < 0.10, "and is already within ten per cent further out")
  note(o$X * Afun(o$r) < 0, "X and the anisotropy have OPPOSITE signs")
}
cat("   So a negative anisotropy gives a positive X, and the radial null contraction\n")
cat("   |f|(k^t)^2 X carries X's sign.\n")

cat("\n=== 3b. and the coupling decides which term leads, which changes the statement ===\n")
cat("   The trace is (6 xi - 1)(1/2) box <phi^2> - m^2 <phi^2>. The m^2 piece carries D^{-3/2},\n")
cat("   subleading to the anisotropy's D^{-5/2}; the box piece carries two more derivatives and so\n")
cat("   D^{-7/2}, which is MORE singular than the anisotropy. So which term leads in X = T - 2A - 4Y\n")
cat("   depends on the coupling, and the powers put the ordering on the page:\n")
pw <- c("the anisotropy A" = -2.5, "trace, m^2 piece" = -1.5, "trace, box piece" = -3.5)
for (nm in names(pw)) cat(sprintf("      %-20s D^{%+.1f}\n", nm, pw[[nm]]))
note(pw[["trace, box piece"]] < pw[["the anisotropy A"]] &&
     pw[["the anisotropy A"]] < pw[["trace, m^2 piece"]],
     "the box term is the most singular and the m^2 term the least")
cat("   At xi = 1/6 the box term is absent and section 3 applies: X -> -2A, opposite signs. Away\n")
cat("   from it the box term leads and X follows the trace instead, checked by integrating the same\n")
cat("   system with (6 xi - 1) = 1:\n")
Tbox <- function(r) 0.5*(M - r)^(-3.5)
stepb <- function(r0, Y0, rend, n = 400000) {
  h <- (rend - r0)/n; r <- r0; Y <- Y0
  for (i in 1:n) {
    X <- Tbox(r) - 2*Afun(r) - 4*Y
    Y <- Y + h*((fpofr(r)/(2*fofr(r)))*X + 2*Afun(r)/r)
    r <- r + h
  }
  list(r = r, Y = Y, X = Tbox(r) - 2*Afun(r) - 4*Y)
}
cat("      D = M - r         X            trace          X/trace\n")
for (Dend in c(1e-3, 1e-4, 1e-5)) {
  o <- stepb(0.5, 0, M - Dend)
  cat(sprintf("   %12.1e %14.4e %14.4e %12.6f\n", Dend, o$X, Tbox(o$r), o$X/Tbox(o$r)))
  note(abs(o$X/Tbox(o$r) - 1) < 0.05, "away from conformal coupling X follows the trace")
}
cat("   So the relation between the two contractions is not one fixed number. At xi = 1/6 they carry\n")
cat("   opposite signs with ratio -2; away from it the radial one follows the trace and the contact\n")
cat("   one does not reach the leading behaviour at all.\n")

cat("\n=== 4. what that means, stated without softening it ===\n")
cat("   The computed negative null-null value belongs to the contact direction, which has no t\n")
cat("   component, so what is negative is the ANISOTROPY. At conformal coupling conservation and the\n")
cat("   trace then make X positive, so T_kk along Penrose's own ingoing radial congruence is POSITIVE\n")
cat("   and the null convergence condition is SATISFIED there. Away from conformal coupling X follows\n")
cat("   the trace and the contact value does not reach the leading behaviour at all. Either way the\n")
cat("   reading that the fold removes Penrose's hypothesis does not follow from the computed sign.\n")
cat("   One thing has to be checked before it is believed, and it is the antecedent. The sign is\n")
cat("   calibrated on the Einstein static universe, where the connecting curve from a point to its\n")
cat("   antipode is separated in TIME as well as angle, so its tangent has a nonzero t component\n")
cat("   and the quantity computed there is a mixture of X and the anisotropy rather than the\n")
cat("   anisotropy alone. Whether the transfer to a hole's contact geodesic carries the whole\n")
cat("   quantity or only one term of it is what decides which way this goes, and this file does not\n")
cat("   settle it.\n")
cat("   What the manuscripts should therefore carry, and now do, is the magnitude and not the sign:\n")
cat("   a thickness rests on the coefficient and the power, and neither of those cares which\n")
cat("   combination of the stress the sign belongs to.\n")
if (fail == 0) cat("\n   all checks passed\n") else {
  cat(sprintf("\n   %d CHECK(S) FAILED\n", fail)); quit(status = 1)
}

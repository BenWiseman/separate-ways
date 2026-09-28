#!/usr/bin/env Rscript
# The band is not decoration: it throws out the most natural rival occupation.
#
# WHY. Section 2.3's mass is a CEILING, because the half-angle state is the least-occupied member
# of the Theta-invariant band and M1 falls as I^{-2/5}. A reader's next question is what the band
# actually forbids, and the honest-sounding answer "everything above the floor" is no answer at
# all, since the floor is where the ceiling already sits. So take the one rival a physicist would
# reach for first. The fold's own map at a horizon produces a THERMAL state, a thermofield double
# at tanh r = e^{-beta omega/2}, which for a fermion is Fermi-Dirac. If the same shape were the
# state at the bang, the occupation would be n = 1/(e^{x/tau} + 1) and not the half-angle form.
#
# Both agree at x = 0, where Theta-invariance pins n = 1/2 and the band has zero width, so the
# comparison is not settled in the infrared where everything agrees. It is settled at x of order 1.
#
# WHAT COMES OUT. Matching the second moment of the pair coherence, which is the matching the
# manuscripts already mention as the reason an exponential converges too, forces tau = 0.3899, and
# at that scale the thermal occupation falls BELOW the band's floor around x = 1. So a thermal
# state at the adopted momentum scale is excluded by Theta-invariance, not by taste. The thermal
# family re-enters the band only for tau >= tau*, computed below, and every admissible member of
# it returns a mass under the ceiling.

TOL <- 1e-9
fail <- 0
note <- function(ok, what) { if (!ok) { cat("   *** FAILED:", what, "\n"); fail <<- fail + 1 } }

nmin <- function(x) { u <- exp(-x^2); u/(2*(1 + sqrt(1 - u))) }   # Boyle-Finn-Turok half-angle
                                                                  # state, written without the
                                                                  # cancellation past x = 6
nmax <- function(x) 1 - nmin(x)                          # unitarity's reflection of it
nFD  <- function(x, tau) 1/(exp(x/tau) + 1)
Ifun <- function(n) integrate(function(x) x^2*n(x), 0, Inf,
                              rel.tol = 1e-12, subdivisions = 4000L)$value/pi^2
M1   <- function(I) 491.6 * (I/I0)^(-2/5)                # sec 2.3: M1 propto I^{-2/5}
   # 491.6 is section 2.3's endpoint, taken as the NORMALISATION of this ratio and not
   # computed here, so no claim row may cite it against this script. What is computed is
   # the ratio I/I0 and everything that follows from it.

cat("=== 1. the quadrature the manuscripts quote, reproduced ===\n")
I0 <- Ifun(nmin)
cat(sprintf("   I(half-angle) = %.13f   residual against the manuscript's figure %.1e\n",
            I0, abs(I0 - 0.0127596673634)))   # the figure itself stays out of the output:
                                              # a provenance gate greps this text, and a script
                                              # that prints the paper's digits can certify them
                                              # by quoting them back (benlm/lessons/LESSONS.md)
note(abs(I0 - 0.0127596673634) < 1e-11, "the half-angle production integral reproduces")
cat(sprintf("   and it returns M1 = %.1f PeV, which is the ceiling\n", M1(I0)))
note(abs(M1(I0) - 491.6) < 1e-9, "the normalisation is the ceiling by construction")

cat("\n=== 2. the second moment of the pair coherence, and the tau it demands ===\n")
cat("   The coherence is 4n(1-n): for the half-angle state it is exp(-x^2) exactly, and for\n")
cat("   Fermi-Dirac it is sech^2(x/2tau). Their second moments are 1/2 and pi^2 tau^2/3.\n")
m2 <- function(w) integrate(function(x) x^2*w(x), 0, Inf, rel.tol = 1e-12)$value /
                  integrate(w, 0, Inf, rel.tol = 1e-12)$value
m2_half <- m2(function(x) exp(-x^2))
tau_m2  <- sqrt(3/(2*pi^2))
m2_FD   <- m2(function(x) 1/cosh(x/(2*tau_m2))^2)
cat(sprintf("   half-angle  <x^2> = %.10f   (analytic 1/2)\n", m2_half))
cat(sprintf("   tau from matching  = %.10f   (analytic sqrt(3)/(sqrt2 pi))\n", tau_m2))
cat(sprintf("   Fermi-Dirac <x^2> at that tau = %.10f\n", m2_FD))
note(abs(m2_half - 0.5) < 1e-10, "the Gaussian coherence has second moment one half")
note(abs(m2_FD - 0.5) < 1e-9, "and the matched thermal one has the same")

cat("\n=== 3. and at that tau the thermal state is under the floor ===\n")
cat("        x        n_FD          floor n_min      n_FD - floor\n")
worst <- Inf
for (x in c(0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 2.0)) {
  d <- nFD(x, tau_m2) - nmin(x); worst <- min(worst, d)
  cat(sprintf("   %6.2f %12.6f %16.6f %16.6f %s\n", x, nFD(x, tau_m2), nmin(x), d,
              ifelse(d < 0, "  <-- OUTSIDE THE BAND", "")))
}
XG  <- c(exp(seq(log(1e-5), log(20), length.out = 40000)))
gap <- function(tau) min(nFD(XG, tau) - nmin(XG))
cat(sprintf("   worst excursion over a fine scan: %.6f\n", gap(tau_m2)))
note(gap(tau_m2) < -1e-4, "the second-moment-matched thermal state leaves the band")
cat("   So the state the fold produces at a HORIZON is not the state at the bang, at the bang's\n")
cat("   own momentum scale. Theta-invariance is what says so.\n")

cat("\n=== 4. where the thermal family re-enters, and what mass it returns there ===\n")
tstar <- uniroot(function(t) gap(t), c(0.30, 3), tol = 1e-14)$root
cat(sprintf("   tau*                       = %.9f\n", tstar))
cat(sprintf("   its second moment          = %.6f   against the adopted 0.5\n", pi^2*tstar^2/3))
Istar <- Ifun(function(x) nFD(x, tstar))
Ianal <- tstar^3 * 1.5 * 1.2020569031595943 / pi^2     # tau^3 (3/2) zeta(3) / pi^2
cat(sprintf("   I(tau*) quadrature         = %.10f\n", Istar))
cat(sprintf("   I(tau*) closed form        = %.10f   (tau^3 (3/2) zeta(3) / pi^2)\n", Ianal))
note(abs(Istar - Ianal) < 1e-10, "the thermal production integral has a closed form and matches")
cat(sprintf("   M1 at tau*                 = %.2f PeV,  %.1f per cent under the ceiling\n",
            M1(Istar), 100*(1 - M1(Istar)/491.6)))
note(M1(Istar) < 491.6, "every admissible thermal member sits under the ceiling")
cat("        tau     in the band?        I          M1 (PeV)\n")
for (t in c(0.39, 0.45, tstar, 0.60, 0.80, 1.20)) {
  Ii <- Ifun(function(x) nFD(x, t))
  cat(sprintf("   %8.4f %12s %14.6f %13.1f\n", t,
              ifelse(gap(t) >= -1e-9, "yes", "no"), Ii, M1(Ii)))
  note((t >= tstar - 1e-6) == (gap(t) >= -1e-9), sprintf("tau* is the boundary (tau = %g)", t))
}
cat("   The family has no upper end, so it puts no floor under the mass, and none is claimed.\n")

cat("\n=== 4b. tau* is exactly one half, and the infrared is why ===\n")
cat("   Expand both about x = 0. The band closes there, so what decides admissibility is the\n")
cat("   SLOPE with which each leaves one half:\n")
cat("      floor        n_min = 1/2 - x/2      + x^3/8      + O(x^5)\n")
cat("      Fermi-Dirac  n_FD  = 1/2 - x/(4tau) + x^3/(48 tau^3) + O(x^5)\n")
cat("   so the linear terms agree at tau = 1/2 and at that value the cubic terms leave\n")
cat("   1/(48 tau^3) - 1/8 = 1/6 - 1/8 = 1/24 > 0. The thermal state therefore sits just above\n")
cat("   the floor for small x precisely when tau >= 1/2, and the numerical root confirms it.\n")
xs <- c(0.1, 0.05, 0.025)
e1 <- sapply(xs, function(x) (nmin(x)    - (1/2 - x/2 + x^3/8))/x^5)
e2 <- sapply(xs, function(x) (nFD(x,0.5) - (1/2 - x/2 + x^3/6))/x^5)
cat(sprintf("   floor residual / x^5 at x = .1, .05, .025: %8.5f %8.5f %8.5f   -> -5/192 = %.5f\n",
            e1[1], e1[2], e1[3], -5/192))
cat(sprintf("   FD    residual / x^5 at tau = 1/2, same x: %8.5f %8.5f %8.5f   -> -1/15  = %.5f\n",
            e2[1], e2[2], e2[3], -1/15))
note(abs(e1[3] + 5/192) < 2e-4 && abs(e2[3] + 1/15) < 2e-4,
     "both series are right through the x^3 term, with the x^5 coefficients they should have")
cat(sprintf("   tau* - 1/2 = %.3e, which is the root finder's tolerance and not a shift\n", tstar - 0.5))
note(abs(tstar - 0.5) < 1e-6, "tau* is one half")
cat(sprintf("   at tau = 1/2 the minimum of n_FD - n_min over x > 0 is %.3e, reached as x -> 0\n",
            min(nFD(XG, 0.5) - nmin(XG))))
note(min(nFD(XG, 0.5) - nmin(XG)) > -1e-10, "and nothing else binds at tau = 1/2")
cat("   The plant: tau = 0.499 must fail, and it must fail in the infrared and not elsewhere.\n")
d499 <- nFD(XG, 0.499) - nmin(XG)
cat(sprintf("      worst excursion %.3e at x = %.4f\n", min(d499), XG[which.min(d499)]))
note(min(d499) < 0 && XG[which.min(d499)] < 0.2, "plant: just below one half it fails, and in the infrared")
cat("   So the closed form is exact: I(tau*) = 3 zeta(3)/(16 pi^2), and M1 = 491.6 (I0/I*)^{2/5}.\n")
Iexact <- 3*1.2020569031595943/(16*pi^2)
cat(sprintf("      3 zeta(3)/(16 pi^2) = %.10f against the quadrature %.10f\n", Iexact, Istar))
note(abs(Iexact - Istar) < 5e-8, "the thermal ceiling's production integral is 3 zeta(3)/(16 pi^2)")

cat("\n=== 5. the plants ===\n")
cstar <- min(nFD(XG, tau_m2)/nmin(XG))
cat(sprintf("   the floor would have to be scaled to %.4f of itself, at x = %.3f, for the matched\n",
            cstar, XG[which.min(nFD(XG, tau_m2)/nmin(XG))]))
cat("   thermal state to be admitted, so the exclusion is by 31 per cent and not by a hair\n")
note(min(nFD(XG, tau_m2) - 0.69*nmin(XG)) >= 0 && min(nFD(XG, tau_m2) - 0.70*nmin(XG)) < 0,
     "plant: a lowered floor lets the excluded state back in, and the crossing is where it says")
cat(sprintf("   tau* - 1e-6 outside: %s     tau* + 1e-6 inside: %s\n",
            ifelse(gap(tstar - 1e-6) < 0, "yes", "NO"), ifelse(gap(tstar + 1e-6) >= 0, "yes", "NO")))
note(gap(tstar - 1e-6) < 0 && gap(tstar + 1e-6) >= 0, "plant: tau* is a genuine boundary")
cat(sprintf("   a 1 per cent error in I moves M1 by %.2f PeV, so the digits quoted are real\n",
            abs(M1(Istar*1.01) - M1(Istar))))

cat("\n=== 6. the numbers the manuscripts may quote, scaled for the digit checker ===\n")
cat(sprintf("   tau star, times 1e4:              %.0f\n", tstar*1e4))
cat(sprintf("   M1 at tau star in PeV:            %.1f\n", M1(Istar)))
cat(sprintf("   matched tau, times 1e4:           %.0f\n", tau_m2*1e4))
cat(sprintf("   floor excursion at matched tau, times 1e3: %.1f\n", abs(gap(tau_m2))*1e3))
cat(sprintf("   the thermal family's own neutrino endpoint, PeV: %.1f\n", M1(Istar)/2))
xback <- uniroot(function(x) nFD(x, tau_m2) - nmin(x), c(1.5, 3), tol = 1e-12)$root
cat(sprintf("   matched tau:                      %.4f\n", tau_m2))
cat(sprintf("   back inside the band at x =       %.2f\n", xback))
cat(sprintf("   worst excursion at x =            %.2f\n", XG[which.min(nFD(XG, tau_m2)/nmin(XG))]))
cat(sprintf("   per cent under the floor there:   %.0f\n", 100*(1 - cstar)))
cat(sprintf("   I at tau* = 3 zeta(3)/16 pi^2:    %.5f\n", Istar))
cat(sprintf("   M1 ceiling of the thermal family, PeV: %.1f\n", M1(Istar)))

if (fail == 0) cat("\n   all checks passed\n") else {
  cat(sprintf("\n   %d CHECK(S) FAILED\n", fail)); quit(status = 1)
}

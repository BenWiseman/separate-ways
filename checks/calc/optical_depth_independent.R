#!/usr/bin/env Rscript
# optical_depth_independent.R -- re-derive the mass-independence of the interior
# optical depth FROM SCRATCH, symbolically first and numerically second, because
# it is the companion's best candidate for a headline number and it should not be
# quoted on the strength of one script agreeing with itself.
#
# Claim under test: tau_opt = n sigma c t_infall is independent of the hole's mass.

G <- 6.67430e-11; cc <- 2.99792458e8; Msun <- 1.98892e30
mp <- 1.67262192e-27; sig <- 1e-29; eta <- 0.1

cat("=== 1. the scaling argument, checked term by term ===\n")
cat("   Each factor's power of M is measured numerically rather than asserted,\n")
cat("   by evaluating at two masses a decade apart and taking the log slope.\n\n")
slope <- function(f) { a <- f(1e3); b <- f(1e4); log10(b/a)/1 }
t_inf <- function(Ms) pi*G*(Ms*Msun)/cc^3
r_s   <- function(Ms) 2*G*(Ms*Msun)/cc^2
V_int <- function(Ms) r_s(Ms)^3
Mdot  <- function(Ms) (1.26e31*Ms)/(eta*cc^2)
n_int <- function(Ms) (2*Mdot(Ms)*t_inf(Ms)/V_int(Ms))/mp
rate  <- function(Ms) n_int(Ms)*sig*cc
tau_o <- function(Ms) rate(Ms)*t_inf(Ms)

for (nm in c("t_infall","r_s","V_interior","Mdot","n","rate","optical depth")) {
  f <- switch(nm, "t_infall"=t_inf, "r_s"=r_s, "V_interior"=V_int,
              "Mdot"=Mdot, "n"=n_int, "rate"=rate, tau_o)
  cat(sprintf("   %-14s power of M = %+6.3f\n", nm, slope(f)))
}

cat("\n=== 2. validation: the check must be able to fail ===\n")
s <- slope(tau_o)
cat(sprintf("   optical depth slope is %.3e, must be zero for mass-independence\n", s))
stopifnot(abs(s) < 1e-9)
# plant a failure: break the volume scaling and confirm the slope moves off zero
V_bad <- function(Ms) r_s(Ms)^2.5
n_bad <- function(Ms) (2*Mdot(Ms)*t_inf(Ms)/V_bad(Ms))/mp
tau_bad <- function(Ms) n_bad(Ms)*sig*cc*t_inf(Ms)
cat(sprintf("   with the interior volume scaling broken to r_s^2.5: slope = %+.3f (must be nonzero)\n",
            slope(tau_bad)))
stopifnot(abs(slope(tau_bad)) > 0.1)
cat("   the check can say no\n")

cat("\n=== 3. the value, across the astrophysical range ===\n")
cat(sprintf("   %12s %16s %16s\n", "M (Msun)", "optical depth", "annihilates?"))
for (Ms in c(24, 1e3, 1e5, 1e7, 1e9, 1e11)) {
  td <- tau_o(Ms)
  cat(sprintf("   %12.0e %16.2f %16s\n", Ms, td, if (td > 1) "yes" else "no"))
}

cat("\n=== 4. how robust is it to the two inputs that are not geometry? ===\n")
cat("   The cancellation is geometric and cannot be tuned away. The VALUE depends on\n")
cat("   the cross-section and the accretion efficiency, so vary both by a decade:\n\n")
cat(sprintf("   %14s", "sigma (m^2)"))
for (e in c(0.03, 0.1, 0.3)) cat(sprintf(" %12s", sprintf("eta=%.2f", e)))
cat("\n")
for (sg in c(1e-30, 1e-29, 1e-28)) {
  cat(sprintf("   %14.0e", sg))
  for (e in c(0.03, 0.1, 0.3)) {
    sig <<- sg; eta <<- e
    cat(sprintf(" %12.1f", tau_o(1e7)))
  }
  cat("\n")
}
sig <- 1e-29; eta <- 0.1

cat(sprintf("
=== 5. flatly ===

  The cancellation is exact and it is geometry, not coincidence: Eddington
  accretion gives Mdot proportional to M, the infall time from horizon to
  singularity goes as M, and the interior volume as M cubed, so the density falls
  as 1/M, the reaction rate falls as 1/M, and the rate times the time is flat. The
  measured slope is %.1e, and breaking the volume scaling on purpose moves it to
  %+.2f, so the check is not blind.

  The value, %.1f at every mass, is not equally robust. It rides on an
  annihilation cross-section and an accretion efficiency, and across a decade
  either way it runs from single digits to the high hundreds. What survives that
  range is the conclusion rather than the number: the optical depth exceeds unity
  everywhere in the box, so the streams annihilate before the singularity at any
  mass and any plausible cross-section, and the interior of such a hole is
  radiation rather than baryons.

  So the quotable statement is the independence, not the 46.7. A number that does
  not care whether the hole is a stellar remnant or a quasar is worth saying; a
  number that moves by two orders under its own inputs is not worth quoting to
  three figures.
", s, slope(tau_bad), tau_o(1e7)))

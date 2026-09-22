# F-AF found the endpoint falsifier is one event wide, and that the weak link is the
# ASSIGNMENT of an event to the decay component rather than its energy. 5.2 leaves
# flavour free, so flavour cannot discriminate. The remaining handle is DIRECTION: a
# decaying halo traces the column density of dark matter, an astrophysical population
# does not.
#
# Decay flux is proportional to the D-factor, the line-of-sight integral of rho (NOT
# rho^2, which is annihilation), so the anisotropy is milder than the annihilation case
# and that is the point to check rather than assume.
#
# NFW: rho(r) = rho_s / [(r/r_s)(1 + r/r_s)^2], r_s = 20 kpc, anchored to
# rho(r_sun) = 0.4 GeV/cm^3 at r_sun = 8.5 kpc.

rs <- 20; rsun <- 8.5; rho_loc <- 0.4
nfw_shape <- function(r) 1/((r/rs)*(1+r/rs)^2)
rho_s <- rho_loc / nfw_shape(rsun)
rho <- function(r) rho_s * nfw_shape(pmax(r,1e-3))
Dfac <- function(psi, lmax=200) {
  integrate(function(l) rho(sqrt(rsun^2 + l^2 - 2*rsun*l*cos(psi))),
            0, lmax, subdivisions=4000, rel.tol=1e-8)$value }

cat("=== 1. the halo anisotropy a decay signal actually carries\n\n")
cat("        angle from GC     D-factor (GeV/cm^3 kpc)     relative to anti-centre\n")
D180 <- Dfac(pi)
for (deg in c(0, 10, 30, 60, 90, 120, 180)) {
  D <- Dfac(deg*pi/180)
  cat(sprintf("   %14.0f %27.2f %26.2f\n", deg, D, D/D180))
}
cat(sprintf("\n  GC / anti-centre at exactly zero degrees = %.1f, but that number is a CUSP\n", Dfac(0)/D180))
cat("  ARTEFACT: the NFW density goes as 1/r at small r, so a line of sight through the\n")
cat("  centre picks up a logarithmically divergent piece and the value depends on the\n")
cat("  inner cutoff imposed. It should not be quoted as the contrast. The robust\n")
cat("  numbers are the solid-angle averages in part 2, which do not pass through the\n")
cat("  cusp with zero impact parameter. A first draft of this line called 38 'a factor\n")
cat("  of a few', which is wrong twice over.\n")

cat("\n=== 2. so how many events to tell halo decay from isotropy?\n\n")
cat("   Treat the sky as two hemispheres, toward and away from the GC. Under isotropy\n")
cat("   the split is 50/50; under halo decay it is set by the hemisphere-averaged\n")
cat("   D-factors. Binomial test, how many events for a 3-sigma separation.\n\n")
# n was 400, which is NOT converged: the solid-angle mean still moves in the fourth decimal,
# and the 3-sigma count sits exactly on an integer boundary there, giving 58 at n=400 and 57
# from n=800 upward. extragalactic_dilution.R used n=800 and silently disagreed with this
# file by one event. p_near converges to 0.698822 by n=12800; 3200 reproduces every printed
# digit below.
hemi <- function(a, b, n=3200) {                 # solid-angle weighted mean D over a cone band
  ps <- seq(a, b, length.out=n)
  sum(sapply(ps, Dfac)*sin(ps))/sum(sin(ps)) }
Dnear <- hemi(0, pi/2); Dfar <- hemi(pi/2, pi)
p_near <- Dnear/(Dnear+Dfar)
cat(sprintf("   hemisphere-mean D toward GC = %.2f, away = %.2f\n", Dnear, Dfar))
cat(sprintf("   so halo decay puts a fraction %.4f of events in the near hemisphere,\n", p_near))
cat(sprintf("   against %.4f for isotropy. The excess to detect is %.4f.\n\n", 0.5, p_near-0.5))
for (k in c(3, 5)) {
  n <- ceiling((k*0.5/(p_near-0.5))^2)
  cat(sprintf("   events needed for %d-sigma hemispheric separation: %s\n", k, format(n, big.mark=",")))
}

cat("\n=== 3. FLATLY, and this is a problem with the falsifier as stated\n\n")
cat("  The single-event test of 5.1 fires only if the event is SECURELY ASSIGNED to the\n")
cat("  decay component. Flavour cannot assign it, since 5.2 leaves flavour free.\n")
cat("  Direction can, but only statistically, and the number is 57 events for three\n")
cat("  sigma and 159 for five. A first draft of this paragraph said 'order a thousand',\n")
cat("  contradicting the table directly above it by a factor of seventeen.\n")
cat("\n  With the right number the conclusion changes and improves. 57 events is NOT out\n")
cat("  of reach: it is the same order as the 100 to 400 that F-AF found the ENDPOINT\n")
cat("  population test needs. So the two requirements are comparable rather than\n")
cat("  circular, and both are met by one population of order 10^2 to 10^3 events at\n")
cat("  these energies. That is a coherent statement of what the falsifier costs, and it\n")
cat("  is a next-generation-detector number rather than an impossible one.\n")

cat("\n=== 4. the route, same turn\n\n")
cat("  One discriminant survives and is not statistical: COINCIDENCE. An astrophysical\n")
cat("  event can be associated with a source - a blazar flare, a tidal disruption - by\n")
cat("  time and direction together, and the decay component cannot be, because a halo\n")
cat("  has no transients. So the assignment that the falsifier needs is not 'this event\n")
cat("  is from the halo' but the weaker and reachable 'this event has NO astrophysical\n")
cat("  counterpart'. That is a statement multimessenger follow-up already produces per\n")
cat("  event, and it is the form in which the falsifier can actually be run.\n")

# ---------------------------------------------------------------------------------------
# 3.2 also quotes an EXACT binomial pair against the normal approximation above. Nothing in
# the repo computed it until now, so it is computed here: the smallest n at which the count
# expected under p_near is significant at the k-sigma level under p = 1/2, one-sided.
exact_n <- function(p, k) {
  a <- pnorm(-k); n <- 2
  while (n <= 5000 && pbinom(round(n*p) - 1, n, 0.5, lower.tail = FALSE) > a) n <- n + 1
  n }
cat("\n=== 4. the exact binomial, against the normal approximation\n\n")
cat("      k     normal    exact\n")
for (k in c(3, 5))
  cat(sprintf("   %4d %9d %8d\n", k, ceiling((k*0.5/(p_near-0.5))^2), exact_n(p_near, k)))
cat("\n  The normal approximation lands within one event at three sigma and two at five,\n")
cat("  which is adequate at this precision.\n")

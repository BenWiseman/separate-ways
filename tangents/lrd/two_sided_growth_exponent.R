# Section 4.3 ends "We have no test that separates the classes and do not propose one."
# That is a kill with no route out of it, so look for the route.
#
# The fold splits black holes into one-sided (formed by collapse, no past singularity)
# and two-sided (a past singularity, two exteriors). Section 4.3 already establishes the
# accretion statement: gravity is set by the SHARED total mass, radiation in one exterior
# cannot reach gas in the other, so each side may accrete at the full Eddington limit for
# the total, and the growth enhancement is kappa = 1 + Mdot_other/Mdot_ours.
#
# What nobody has done is notice WHERE kappa sits. Eddington-limited growth is
#     dM/dt = kappa * M / t_S      =>     M(t) = M_seed exp(kappa t / t_S),
# because the Eddington rate on each side is set by the total mass. kappa is in the
# EXPONENT, not a prefactor. Over the ~13 e-folds the little red dots need, a few per
# cent of extra accretion from the far side is an order of magnitude in final mass.
#
# That cuts both ways and both are worth having: it is a mechanism for overmassive
# early holes, and it is so sensitive that the observed overmassive factor MEASURES
# kappa rather than merely tolerating it.

H0 <- 67.4; Om <- 0.315; OL <- 0.685
Mpc <- 3.0857e22; Gyr <- 3.1557e16
tz <- function(z) (2/(3*H0*sqrt(OL)))*asinh(sqrt(OL/Om)*(1+z)^(-1.5)) * (Mpc/1e3) / Gyr
t_S <- 45          # Myr, the Salpeter time at eps=0.1, lambda_Edd=1; lrd_budget.R's convention

efolds <- function(zs, zo, kappa=1) (tz(zo)-tz(zs))*1000/t_S*kappa

cat("=== 0. validation: reproduce the number already in the paper\n\n")
N1 <- efolds(20, 7)
cat(sprintf("   z=20 to z=7 at kappa=1 : %.2f e-folds\n", N1))
cat(sprintf("   section 3.4 quotes      : 12.9\n"))
cat(sprintf("   agree to               : %.3f\n", abs(N1-12.9)))
stopifnot(abs(N1 - 12.9) < 0.05)
cat("   -> the time base and Salpeter convention match the paper. Proceed.\n")

cat("\n=== 1. kappa is in the exponent, so price it\n\n")
cat("   Seed mass required to reach 1e7 Msun by z=7, starting at z=20.\n\n")
cat("      kappa   e-folds    seed needed (Msun)   vs kappa=1\n")
base <- 1e7/exp(efolds(20,7,1))
for (k in c(1.0, 1.1, 1.2, 1.5, 2.0)) {
  N <- efolds(20, 7, k); s <- 1e7/exp(N)
  cat(sprintf("   %7.2f %9.2f %20.3g %12.3g\n", k, N, s, s/base))
}
cat("\n   A ten per cent contribution from the far side cuts the required seed by a\n")
cat("   factor of 3.6. Doubling the accretion removes the seed problem outright.\n")

cat("\n=== 2. the same exponent run forwards: how overmassive does it make a hole?\n\n")
cat("   At fixed seed and fixed time, a two-sided hole outweighs a one-sided one by\n")
cat("   exp((kappa-1) * N_1), with N_1 the one-sided e-fold count over the window.\n\n")
cat("      window        N_1     kappa=1.1  kappa=1.2  kappa=1.5   kappa=2\n")
for (w in list(c(20,7), c(20,5), c(15,7), c(10,7))) {
  N <- efolds(w[1], w[2])
  cat(sprintf("   z=%2.0f to z=%1.0f %8.2f", w[1], w[2], N))
  for (k in c(1.1,1.2,1.5,2.0)) cat(sprintf(" %10.3g", exp((k-1)*N)))
  cat("\n")
}

cat("\n=== 3. invert it: the OBSERVED overmassive factor measures kappa\n\n")
cat("   Little red dots sit high on the M_BH-M_star relation. Taking the reported\n")
cat("   excess as a range rather than a number, and attributing ALL of it to the far\n")
cat("   side (an upper bound on kappa, since other mechanisms also lift the ratio):\n\n")
N <- efolds(20, 7)
cat("      excess over local relation    implied kappa    implied Mdot_other/Mdot_ours\n")
for (R in c(3, 10, 30, 100)) {
  k <- 1 + log(R)/N
  cat(sprintf("   %22.0fx %16.3f %28.3f\n", R, k, k-1))
}
cat("\n   Read at face value the whole observed range corresponds to the far side\n")
cat("   accreting at 9 to 36 per cent of our rate. Section 3b shows that reading is\n")
cat("   premature: it assumes continuous Eddington accretion.\n")

cat("\n=== 3b. the confound that nearly sank section 3: the duty cycle\n\n")
cat("   None of the above survives if holes accrete at Eddington only part of the\n")
cat("   time. With duty cycle lambda the one-sided e-fold count is lambda*N_1, and\n")
cat("   the overmassive factor is exp((kappa-1) lambda N_1). So what the excess\n")
cat("   measures is the PRODUCT (kappa-1) lambda, not kappa. Price the degeneracy:\n\n")
cat("      lambda   lambda*N_1    kappa needed for 10x   kappa needed for 100x   is kappa=2 excluded?\n")
for (lam in c(1.0, 0.5, 0.3, 0.1, 0.03)) {
  Ne <- lam*N
  k10 <- 1 + log(10)/Ne; k100 <- 1 + log(100)/Ne
  excl <- if (exp(1.0*Ne) > 100) "yes" else "NO"
  cat(sprintf("   %8.2f %12.2f %22.2f %23.2f %22s\n", lam, Ne, k10, k100, excl))
}
cat("\n   At lambda=1 a far side matching ours gives 4e5, which is excluded. At\n")
cat("   lambda=0.1 the same kappa=2 gives only a factor of 3.6, which is squarely\n")
cat("   inside what is observed. So the exclusion holds ONLY for near-continuous\n")
cat("   Eddington accretion, and the literature pressure on these objects is\n")
cat("   precisely that they cannot sustain that.\n")

cat("\n=== 4. the discriminator section 4.3 says it does not have\n\n")
cat("   Soltan's argument compares the local black hole mass density with the mass\n")
cat("   accreted, inferred from integrated AGN light at radiative efficiency eta:\n")
cat("       rho_accreted = (1-eta)/(eta c^2) * integral of L dt.\n")
cat("   Mass that arrived from the far side emitted its light into the far side, so it\n")
cat("   is in rho_BH and not in the integral. A survey fitting eta to a population with\n")
cat("   a two-sided fraction f_2s recovers\n")
cat("       eta_apparent = eta_true / (1 + f_2s (kappa-1)),\n")
cat("   i.e. it UNDERSTATES the efficiency. Run it:\n\n")
eta_true <- 0.10
cat("      f_2s   kappa    eta_apparent    deficit\n")
for (f in c(0.1, 0.3, 1.0)) for (k in c(1.2, 1.36)) {
  ea <- eta_true/(1+f*(k-1))
  cat(sprintf("   %6.2f %7.2f %15.4f %9.1f%%\n", f, k, ea, 100*(ea/eta_true-1)))
}
cat("\n   Published Soltan efficiencies sit near 0.07-0.10 with systematics of tens of\n")
cat("   per cent from the bolometric correction and the obscured fraction alone. A\n")
cat("   4-26 per cent deficit is INSIDE that, so the test does not fire today.\n")

cat("\n=== 5. flatly\n\n")
cat("  FAILS, twice, and the second one is mine:\n\n")
cat("   (i)  This does not separate the classes with present data. The Soltan deficit\n")
cat("        a plausible kappa produces is smaller than the bolometric-correction\n")
cat("        systematic, so section 4.3's sentence stands as written.\n\n")
cat("   (ii) The overmassive factor does NOT measure kappa. It measures the product\n")
cat("        (kappa-1)*lambda*N_1, and the duty cycle lambda is not known for these\n")
cat("        objects. Section 3 of this script claimed a measurement; section 3b\n")
cat("        withdraws it. kappa=2 is excluded only for lambda above about 0.3, and\n")
cat("        the literature's whole complaint about these sources is that they cannot\n")
cat("        sustain lambda=1. Worse, the degeneracy runs the wrong way for the fold:\n")
cat("        at low duty cycle ANY excess can be fitted by raising kappa, so the\n")
cat("        mechanism becomes flexible rather than predictive. A model that can fit\n")
cat("        any overmassive ratio has not explained the overmassive ratio.\n\n")
cat("  WHAT SURVIVES:\n")
cat("   (a) The structural point is untouched and is not in the paper. Section 4.3\n")
cat("       treats the far side as a factor of two on the RATE. It is a factor in the\n")
cat("       EXPONENT of the growth law, because both sides are Eddington-limited by\n")
cat("       the shared mass. Over the little red dots' window that converts a factor\n")
cat(sprintf("       of two into %.2g at full duty cycle. The paper understates its own\n", exp(N)))
cat("       mechanism by five orders of magnitude.\n")
cat("   (b) The fold needs no new parameter to do this. kappa is fixed by the far\n")
cat("       side's gas supply, which is not ours to choose, and one-sided holes have\n")
cat("       kappa=1 identically.\n\n")
cat("  NEXT ROUTE, and it is degeneracy-free, which is why it is worth stating:\n")
cat("  kappa is not a continuous knob across the population. Section 4.3 argues that\n")
cat("  only non-collapse seeds can be two-sided, and a collapse-formed hole has\n")
cat("  kappa=1 EXACTLY, not approximately. Duty cycle, by contrast, is continuous and\n")
cat("  shared by both classes. So the fold predicts the M_BH-M_star residual is drawn\n")
cat("  from TWO sequences rather than one broadened distribution, with the split set\n")
cat("  by seed origin and not by accretion history. Bimodality survives the duty-cycle\n")
cat("  degeneracy because lambda smears each sequence without merging them. Testing it\n")
cat("  needs a residual distribution rather than a mean, which is a catalogue\n")
cat("  calculation and a companion paper. A dip test or a two-component mixture fit on\n")
cat("  the published residuals would do it, and that is the cheapest real test this\n")
cat("  framework has yet offered for the two-sided class.\n")

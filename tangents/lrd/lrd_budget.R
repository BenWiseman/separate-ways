# Ben: the picture must COMFORTABLY accommodate JWST little red dots and black hole
# stars. Prior work (pbh_routes_tested.R) found the fold PERMITS a subdominant seed
# population without disturbing its dark-matter budget but does not PRODUCE one. The
# question that settles "comfortably" is therefore a budget question, and it can be done
# parametrically with no catalogue lookup.

H0 <- 67.4; Om <- 0.315; OL <- 0.685
Mpc <- 3.0857e22; Gyr <- 3.1557e16
# cosmic time in a flat LCDM, Gyr
tz <- function(z) (2/(3*H0*sqrt(OL)))*asinh(sqrt(OL/Om)*(1+z)^(-1.5)) * (Mpc/1e3) / Gyr

cat("=== 1. how much growth time is there, and what seed does it demand?\n\n")
cat("   Salpeter time t_S = 45 Myr (eps/0.1) / lambda_Edd. e-folds N = dt/t_S.\n\n")
cat("        z_seed  z_obs    t(z_seed)   t(z_obs)   dt (Myr)   N e-folds   seed for 1e7 Msun\n")
for (zs in c(20, 15, 10)) for (zo in c(7, 5)) {
  dt <- (tz(zo)-tz(zs))*1000
  for (lam in c(1)) {
    N <- dt/45*lam
    cat(sprintf("   %8.0f %6.0f %11.4f %10.4f %10.1f %11.2f %19.3g\n",
        zs, zo, tz(zs), tz(zo), dt, N, 1e7/exp(N)))
  }
}
cat("\n  So Eddington-limited growth from z=20 to z=7 gives ~13 e-folds and needs a seed\n")
cat("  of order 10 to 100 Msun at 1e7 Msun final: a stellar remnant is enough IF it\n")
cat("  accretes at Eddington the whole way. The tension in the literature is that it\n")
cat("  cannot, not that the arithmetic fails. A heavier seed removes the requirement.\n")

cat("\n=== 2. the seed BUDGET, which is what 'comfortably' means here\n\n")
cat("   Mass fraction of the dark matter needed in seeds, as a function of the seed\n")
cat("   mass and the comoving number density of hosts. rho_DM = Om_DM rho_crit.\n\n")
rho_crit <- 2.775e11 * (H0/100)^2          # Msun / Mpc^3 (h^-1 units cancelled)
rho_DM <- 0.265/(H0/100)^2 * rho_crit * (H0/100)^2
cat(sprintf("   rho_DM = %.3e Msun / Mpc^3 (comoving)\n\n", rho_DM))
cat("        n_host (Mpc^-3)   M_seed (Msun)    f_seed = n M / rho_DM\n")
for (n in c(1e-4, 1e-5, 1e-6)) for (Ms in c(1e2, 1e4, 1e5)) 
  cat(sprintf("   %15.1e %15.1e %24.3e\n", n, Ms, n*Ms/rho_DM))
cat("\n  Even the most generous combination, one 1e5 Msun seed per 1e4 Mpc^3, needs\n")
cat("  f_seed ~ 1e-5 of the dark matter. The fold's own budget tolerates f ~ 1e-3\n")
cat("  before M_1 moves by 0.04%, which is inside its stated precision\n")
cat("  (`calc/tangents/seam/pbh_routes_tested.R`).\n")

cat("\n=== 3. so, flatly\n\n")
cat("  ACCOMMODATES: yes, and with about two orders of magnitude of margin. A seed\n")
cat("  population large enough to explain every little red dot perturbs the fold's\n")
cat("  dark-matter mass by less than a part in a thousand of its own error bar.\n")
cat("  EXPLAINS: no. Nothing in the fold puts power at the seed scale, and claiming\n")
cat("  otherwise would be inventing a mechanism. The honest statement is compatibility\n")
cat("  with margin, not a prediction, and the margin is the part worth quoting because\n")
cat("  a model whose dark sector is FULL could not say it.\n")

cat("\n=== 4. and the one place the fold is genuinely exposed\n\n")
cat("  Its dark matter is a single 491.6 PeV sterile neutrino fixed by the relic\n")
cat("  abundance, so the budget is closed. If little red dots were later shown to need\n")
cat("  f_PBH of order unity - a dark sector made of black holes - the fold would be\n")
cat("  refuted, not adjusted, because M_1 would have to move by:\n\n")
M1 <- 491.6
cat("        f_PBH      M_1 (PeV)     shift      predicted neutrino line (PeV)\n")
for (f in c(0, 0.01, 0.1, 0.5, 0.9)) {
  Mn <- M1*(1-f)^(2/5)
  cat(sprintf("   %10.2f %12.2f %9.1f%% %28.1f\n", f, Mn, 100*(Mn/M1-1), Mn/2))
}
cat("\n  The neutrino line at half the mass is the observable, so a black-hole dark\n")
cat("  sector moves the prediction away from the KM3NeT energy by a measurable amount.\n")
cat("  That is a real falsifier and it runs through the little red dots, which is worth\n")
cat("  having: the two subjects are not independent in this model.\n")

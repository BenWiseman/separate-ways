# C.2 asserts two numbers with no computation behind either of them:
#
#   eta_* = 4.39/H0          the conformal time at which a horizon identification would act
#   eta_prod/eta_* = 1.6e-25 production being that much earlier
#
# OPEN_FOR_BEN_20260922.md section 3 filed both as unreproducible and named recomputing them
# as the next action. This is that recomputation. Base R only, per the repo rule.
#
# The route, which the paper states but does not carry out:
#   a(eta) = a_1 eta in radiation domination, so a_1 = H0 sqrt(Omega_r) with a = 1 today
#   gamma  = M_1 a_1                              (2.2's definition)
#   eta_prod ~ 1/sqrt(gamma)                      (C.2's own statement)
#   eta_*  = integral of da/(a^2 E(a)) from 0 to infinity

H0_eV    <- 1.437e-33      # 67.4 km/s/Mpc in eV
M1_eV    <- 4.916e17       # 491.6 PeV, the paper's ceiling
Om       <- 0.315
OL       <- 0.685
Og       <- 5.43e-5        # photons alone, Omega_gamma h^2 = 2.47e-5, h = 0.674
Or_full  <- 9.24e-5        # photons plus three massless neutrino species

# --- eta_*: total conformal time, bang to a -> infinity, in units of 1/H0
eta_total <- function(Or, a_hi = 1e8, n = 400000) {
  f <- function(a) 1/(a^2 * sqrt(Om*a^-3 + Or*a^-4 + OL))
  la <- log(1e-14); lb <- log(a_hi); h <- (lb-la)/n
  i <- 0:n; a <- exp(la + i*h)
  w <- ifelse(i == 0 | i == n, 1, ifelse(i %% 2 == 1, 4, 2))
  sum(w * f(a) * a) * h/3                      # da = a dln a
}

cat("=== 1. eta_*, the total conformal time\n\n")
cat(sprintf("   matter + Lambda only     : %.4f / H0\n", eta_total(0)))
cat(sprintf("   with photons             : %.4f / H0\n", eta_total(Og)))
cat(sprintf("   with photons + neutrinos : %.4f / H0\n", eta_total(Or_full)))
cat("\n   The paper's 4.39 is the matter + Lambda value. Radiation shortens it to 4.33,\n")
cat("   so the quoted figure neglects radiation, which for a conformal time running to\n")
cat("   infinity is a 1.4 per cent choice and not a mistake.\n")

# --- eta_prod/eta_*
cat("\n=== 2. eta_prod / eta_*, production against that boundary\n\n")
cat("        Omega_r used       a_1 [eV]      gamma [eV^2]   eta_prod/eta_*\n")
for (nm in c("photons only", "photons + nu")) {
  Or <- if (nm == "photons only") Og else Or_full
  a1 <- H0_eV * sqrt(Or)                       # a = a_1 eta in radiation domination
  gm <- M1_eV * a1
  ratio <- (1/sqrt(gm)) / (eta_total(Or)/H0_eV)
  cat(sprintf("   %-16s %12.3e %14.3e %14.2e\n", nm, a1, gm, ratio))
}

cat("\n=== 3. flatly\n\n")
cat("  eta_* = 4.39/H0 REPRODUCES exactly, as the matter-plus-Lambda total conformal time.\n")
cat("  eta_prod/eta_* comes out at 1.4e-25 with photons alone and 1.3e-25 with neutrinos\n")
cat("  included, against the 1.6e-25 the paper quotes. Same order, same route, and the\n")
cat("  residual factor of 1.2 is inside the freedom in Omega_r and in reading '~' in\n")
cat("  eta_prod ~ 1/sqrt(gamma). The paper now quotes a rounded 1e-25 and cites this file\n")
cat("  rather than a digit it cannot reproduce.\n\n")
cat("  NOT RESOLVED HERE: the companion assertion that the largest conceivable\n")
cat("  holonomy-sector differential in n_dm/s is 5.8e-26. Nothing in the paper relates it\n")
cat("  to the ratio above, the quotient 1.6/0.58 = 2.76 matches no constant in reach, and\n")
cat("  no derivation is shown. It has been cut from C.2 rather than carried unsourced.\n")

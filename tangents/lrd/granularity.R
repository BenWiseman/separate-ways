# Two questions about the little red dots that the paper has not asked.
# (1) A 491.6 PeV relic is absurdly heavy per particle, so the dark matter is GRANULAR.
#     How granular, in units a person can picture?
# (2) Discrete dark matter seeds structure by Poisson noise, and the effect GROWS with the
#     particle mass. Superheavy dark matter is exactly the case where this is largest. Does
#     it reach the seed masses the little red dots need? Compute, do not hope.

eV_kg <- 1.78266192e-36
Msun  <- 1.98892e30
M1    <- 491.6e15 * eV_kg                    # kg per particle
rho   <- 0.12 * 1.878e-26                    # kg/m^3, Omega_DM h^2 = 0.12
n     <- rho/M1
cat("=== 1. how granular is this dark matter?\n\n")
cat(sprintf("   mass per particle      %.3e kg  (a grain of sand is ~1e-8 kg)\n", M1))
cat(sprintf("   number density         %.3e per m^3\n", n))
cat(sprintf("   mean spacing           %.1f m\n", n^(-1/3)))
cat(sprintf("   particles per km^3     %.2f\n", n*1e9))
cat(sprintf("   particles in the Earth's volume (1.08e21 m^3): %.3e\n", n*1.08e21))
cat("\n   For comparison, a 100 GeV WIMP at the same mass density:\n")
mW <- 100e9*eV_kg; nW <- rho/mW
cat(sprintf("   spacing %.2f m, i.e. this model's dark matter is %.0f times sparser in\n",
    nW^(-1/3), (n^(-1/3))/(nW^(-1/3))))
cat("   number and correspondingly heavier per particle.\n")

# --- WHICH density? This mattered and was wrong in the manuscript.
# Everything above uses Omega_DM h^2 = 0.12, which is the COSMOLOGICAL MEAN. 3.4 called it
# "the measured dark-matter density" and then said this picture was behind the event count in
# Appendix D.1 -- but D.1 correctly uses a LOCAL halo density of 0.4 GeV/cm^3, which is
# 3.2e5 times denser. Two different densities were being linked as though they were one.
# Quote both, labelled, and use the local one for any claim about an instrument.
GeV_kg  <- 1.78266192e-27
rho_loc <- 0.4 * GeV_kg * 1e6                # 0.4 GeV/cm^3 -> kg/m^3
n_loc   <- rho_loc/M1
nW_loc  <- rho_loc/mW
cat("\n   Cosmological mean against the local halo, since they differ by 3e5:\n\n")
cat(sprintf("      cosmological mean   %9.2f per km^3, spacing %7.1f m   (WIMP %6.2f m)\n",
    n*1e9, n^(-1/3), nW^(-1/3)))
cat(sprintf("      local halo, 0.4     %9.3g per km^3, spacing %7.1f m   (WIMP %6.3f m)\n",
    n_loc*1e9, n_loc^(-1/3), nW_loc^(-1/3)))
cat(sprintf("      local is %.2e times denser in number.\n", n_loc/n))
cat("\n   The instrument claim belongs to the LOCAL row: even there the particles sit about\n")
cat("   eleven metres apart, against six centimetres for a 100 GeV WIMP, so the statement\n")
cat("   that this dark matter is not a fluid on any scale an instrument spans survives at\n")
cat("   the density a detector actually sees, which is the density D.1 uses.\n")
stopifnot(abs(n*1e9 - 2.57) < 0.05, abs(n^(-1/3) - 730) < 3,
          abs(n_loc*1e9 - 8.14e5) < 1e4, abs(n_loc^(-1/3) - 10.7) < 0.2)

cat("\n=== 2. Poisson seeding: the one mechanism that FAVOURS a heavy particle\n\n")
cat("   A region holding N particles has a fractional shot-noise overdensity 1/sqrt(N), and\n")
cat("   N = M/m, so delta_P = sqrt(m/M). Heavier particle, larger seed noise. Evaluate it.\n\n")
cat("      seed mass M        N particles      delta_Poisson    reaches 1e-3?\n")
for (Ms in c(1e2,1e4,1e5,1e6,1e8)) {
  M <- Ms*Msun; N <- M/M1; d <- sqrt(M1/M)
  cat(sprintf("   %8.0e Msun %18.3e %16.3e %16s\n", Ms, N, d, if (d>1e-3) "yes" else "no"))
}
Mneed <- M1/(1e-3)^2
cat(sprintf("\n   delta_P reaches 1e-3 only below M = %.3e kg = %.3e Msun.\n", Mneed, Mneed/Msun))
cat("   FLAT NEGATIVE, stated with the right exponent (corrected 2026-09-21: the old line said\n")
cat("   thirty-odd, which matched none of the three comparisons available). Against the 1e-3 a\n")
cat("   seed needs, the shortfall is 23.7 orders; against unity it is 26.7; and in SEED MASS the\n")
cat("   gap to 1e4-1e5 Msun is 47.4 orders. The old line conflated the first with the third.\n")
cat("   the little red dots need. Poisson seeding from this particle is irrelevant, and the\n")
cat("   fact that the mechanism favours heavy dark matter does not rescue it.\n")

cat("\n=== 3. so what CAN the paper say, and it is not nothing\n\n")
rho_seed <- function(Ms, nhost) Ms*Msun*nhost / (rho * (3.0857e22)^3)
for (Ms in c(1e4,1e5,1e6)) {
  f <- rho_seed(Ms, 1e-4)
  cat(sprintf("   seeds of %.0e Msun at 1e-4 per Mpc^3 take %.2e of the dark budget\n", Ms, f))
}
cat("\n   The budget is closed, every gram is the sterile neutrino, and heavy seeds still fit\n")
cat("   with orders to spare. That is the honest claim: not that the fold makes the seeds,\n")
cat("   but that closing the dark sector does not forbid them, which is the thing a closed\n")
cat("   sector is most likely to get wrong.\n")

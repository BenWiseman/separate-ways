# ==========================================================================================
# The extragalactic shape of 3.2 is zero-parameter, so the ONLY free quantity in the decay
# component is 1/tau. That means the absolute flux is predicted once tau is named, and the
# paper has so far treated tau purely as an input. Compute the normalisation.
#
# Two independent pieces, and their ratio is a check on the 42.7 per cent computed in
# extragalactic_dilution.R by an entirely different route (absolute fluxes, not column ratios).
# ==========================================================================================
c_cm  <- 2.99792458e10                 # cm/s
Mpc   <- 3.0856775814913673e24         # cm
kpc   <- Mpc/1e3
H0    <- 67.4/(Mpc/1e5)                # s^-1  (67.4 km/s/Mpc)
Om    <- 0.315; OL <- 0.685
Ez    <- function(z) sqrt(Om*(1+z)^3 + OL)
rhoDM <- 0.1200/0.674^2 * 1.05371e-5*0.674^2    # GeV/cm^3, comoving mean
M1    <- 491.6e6                       # GeV
Enu   <- M1/2                          # GeV, the two-body line
Nnu   <- 0.5                           # hard neutrinos per decay, the 1:1:2 tree-level yield

cat(sprintf("  H0 = %.4e s^-1,  rho_DM = %.4e GeV/cm^3,  M1 = %.4e GeV,  E_nu = %.4e GeV\n",
            H0, rhoDM, M1, Enu))

# --- extragalactic: E^2 dPhi/dE = c rho N E / (4 pi M1 tau H(z*)),  z* = E_nu/E - 1 --------
E2dPhi_EG <- function(E, tau) {
  z <- Enu/E - 1
  c_cm*rhoDM*Nnu*E/(4*pi*M1*tau*H0*Ez(z))
}
# --- Galactic line: integrated number flux (1/4pi) N J / (M1 tau), then weight by E --------
# hemisphere-averaged NFW column from extragalactic_dilution.R, in GeV cm^-3 kpc
Jbar <- 7.1508 * kpc                    # GeV/cm^2
EPhi_gal <- function(tau) Enu * Nnu*Jbar/(4*pi*M1*tau)

for (tau in c(1e26, 1e27, 1e28, 1e29, 1e30)) {
  eg <- E2dPhi_EG(Enu, tau); gal <- EPhi_gal(tau)
  cat(sprintf("\n  tau = %.0e s\n", tau))
  cat(sprintf("     extragalactic, at the line  E^2 dPhi/dE = %.3e GeV cm^-2 s^-1 sr^-1\n", eg))
  cat(sprintf("     Galactic line, E x (integrated flux)    = %.3e GeV cm^-2 s^-1 sr^-1\n", gal))
  cat(sprintf("     extragalactic share                     = %.4f\n", eg/(eg+gal)))
}

# --- is the ratio above an independent check on the 0.427? Work it out rather than claim it.
Ikern <- function(z) 1/((1+z)*Ez(z))
Iz    <- integrate(Ikern, 0, 1000, rel.tol=1e-12)$value
tot_EG  <- function(tau) c_cm*rhoDM*Nnu/(4*pi*M1*tau*H0)*Iz     # integrated over all E
tot_gal <- function(tau) Nnu*Jbar/(4*pi*M1*tau)
cat("\n  ------------------------------------------------------------------------------\n")
cat("  Is this an independent check on the 0.427 of extragalactic_dilution.R? No.\n")
cat("  ------------------------------------------------------------------------------\n")
cat(sprintf("  The share printed above, %.4f, compares a DIFFERENTIAL quantity at the line\n",
            E2dPhi_EG(Enu,1e28)/(E2dPhi_EG(Enu,1e28)+EPhi_gal(1e28))))
cat("  with an INTEGRATED one, so it is not the same ratio and its 3 per cent offset from\n")
cat("  0.427 means nothing either way.\n\n")
cat("  Matched properly, integrating the extragalactic flux over all E:\n")
cat(sprintf("     Int dz/((1+z)E(z))            = %.6f\n", Iz))
cat(sprintf("     total extragalactic / total   = %.6f\n", tot_EG(1e28)/(tot_EG(1e28)+tot_gal(1e28))))
cat("  and that IS 0.427. But it is not an independent route: substituting z = E_nu/E - 1\n")
cat("  turns Int dE/(E H(z*)) into (1/H0) Int dz/((1+z)E(z)) identically, so the flux ratio\n")
cat("  reduces to c rho Iz/(H0 Jbar), which is the column ratio. Same calculation in\n")
cat("  different variables. It checks the algebra, not the physics.\n")

cat("\n  What IS new here is the absolute normalisation, which the paper has never stated.\n")

cat("\n  Scaling, so anyone can confront it with an exposure:\n")
cat(sprintf("     E^2 dPhi/dE (extragalactic, at the line) = %.3e x (1e28 s / tau)\n",
            E2dPhi_EG(Enu, 1e28)))
cat("     in GeV cm^-2 s^-1 sr^-1. The shape away from the line carries no further\n")
cat("     freedom, so a limit at ANY single energy converts to a limit on tau with no\n")
cat("     spectral assumption. That is what the zero-parameter shape buys.\n")

cat("\n  For scale only, NOT as a limit: IceCube's measured per-flavour astrophysical flux\n")
cat("  is E^2 Phi = 1.66e-8 GeV cm^-2 s^-1 sr^-1 at 100 TeV (arXiv:2001.09520, gamma=2.53).\n")
cat("  The prediction at the line is the same order at tau ~ 1e28 s, so the model sits in\n")
cat("  the observationally interesting range rather than far outside it. Turning that into\n")
cat("  a bound needs the exposure at 245.8 PeV, which is not computed here.\n")

# How big is the Nariai hole in OUR universe?
#
# The abstract and section 6 call the qualifying hole "the maximal hole". That is a
# statement about SdS, not about a size. This turns it into a length, because a reader
# needs to know whether the class the fold allows is astrophysical or cosmological.
#
# Degenerate SdS: f(r) = 1 - 2M/r - Lr^2/3 with f = f' = 0.
#   f' = 0  =>  M = L r^3 / 3
#   f  = 0  =>  1 - 2Lr^2/3 - Lr^2/3 = 1 - L r^2 = 0  =>  r0 = 1/sqrt(L)
# Pure de Sitter has r_dS = sqrt(3/L), so r0 = r_dS/sqrt(3) exactly.

cat("=== Nariai radius in physical units ===\n\n")

# Planck 2018 TT,TE,EE+lowE+lensing+BAO
H0    <- 67.36            # km/s/Mpc
OmL   <- 0.6847
cc    <- 299792.458       # km/s
Mpc_ly <- 3.261563777e6   # light years per Mpc

hubble_Mpc <- cc / H0                 # c/H0 in Mpc
hubble_Gly <- hubble_Mpc * Mpc_ly / 1e9

# L = 3 (H0/c)^2 OmL, so 1/sqrt(L) = (c/H0) / sqrt(3 OmL)
r0_Gly  <- hubble_Gly / sqrt(3 * OmL)
rdS_Gly <- hubble_Gly / sqrt(OmL)

cat(sprintf("  Hubble radius c/H0        %10.4f Gly\n", hubble_Gly))
cat(sprintf("  de Sitter radius sqrt(3/L)%10.4f Gly\n", rdS_Gly))
cat(sprintf("  Nariai radius    1/sqrt(L)%10.4f Gly\n", r0_Gly))
cat(sprintf("  ratio r0/r_dS             %10.7f   (exact 1/sqrt(3) = %.7f)\n\n",
            r0_Gly/rdS_Gly, 1/sqrt(3)))

stopifnot(abs(r0_Gly/rdS_Gly - 1/sqrt(3)) < 1e-12)

# Mass, for scale. M = L r0^3/3 = r0/3 in geometric units (G=c=1).
Msun_m <- 1476.6250385    # GM_sun/c^2 in metres
ly_m   <- 9.4607304725808e15
M_geom_m <- (r0_Gly * 1e9 * ly_m) / 3
M_sun    <- M_geom_m / Msun_m
cat(sprintf("  mass M = r0/3             %10.4e solar masses\n", M_sun))
cat(sprintf("  rounded to two figures    %.1f x 10^22 solar masses\n", M_sun/1e22))

# Comparisons a reader can picture.
cat(sprintf("\n  compare: observable universe, comoving radius   46.5   Gly\n"))
cat(sprintf("           light-travel distance to the CMB        45.6   Gly (comoving)\n"))
cat(sprintf("           Hubble radius                           %6.2f Gly\n", hubble_Gly))
cat(sprintf("\n  So the Nariai horizon is about %.0f billion light years in radius:\n", r0_Gly))
cat(sprintf("  a factor %.1f short of the comoving radius of the visible universe,\n", 46.5/r0_Gly))
cat(sprintf("  and %.2f times the Hubble radius.\n", r0_Gly/hubble_Gly))
cat(sprintf("\n  CONCLUSION: 'the size of the visible universe' is wrong by a factor of ~%.1f.\n", 46.5/r0_Gly))
cat(sprintf("  The defensible phrase is 'about %.0f billion light years across' for the DIAMETER\n", 2*r0_Gly))
cat(sprintf("  (2 r0 = %.2f Gly), or 'ten billion light years in radius'.\n", 2*r0_Gly))

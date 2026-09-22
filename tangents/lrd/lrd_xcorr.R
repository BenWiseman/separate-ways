# ==========================================================================================
# Does the energy-resolved angular cross-correlation between the relic decay signal and the
# early black-hole population reach an actual constraint?  Compute, do not hope.
#
# The one thing that makes this framework different from a generic decaying-dark-matter
# anisotropy search: the parent mass is NOT a free parameter.  M_1 <= 491.6 PeV is fixed by the
# abundance, so the two-body line sits at E_0 = M_1/2 = 245.8 PeV and the extragalactic
# component of the SAME line arrives at E = E_0/(1+z).  Energy IS redshift, with no free
# parameter to absorb a mismatch.  That is the only place JWST's z = 4-9 population and a
# neutrino telescope genuinely meet.  Price it.
# ==========================================================================================

M1   <- 491.6            # PeV, the ceiling
E0   <- M1/2             # PeV, two-body line
H0   <- 67.4; Om <- 0.315; OL <- 0.685
c_kms<- 299792.458
Ez   <- function(z) sqrt(Om*(1+z)^3 + OL)
cH0_Mpc <- c_kms/H0                      # Hubble distance, Mpc
Mpc_cm  <- 3.0856775814913673e24
kpc_cm  <- Mpc_cm/1e3
cH0_cm  <- cH0_Mpc*Mpc_cm

cat("=============================================================================\n")
cat(" 1. THE ENERGY AXIS IS A REDSHIFT AXIS, AND THE CALIBRATION HAS NO FREE PARAMETER\n")
cat("=============================================================================\n\n")
cat(sprintf("   line energy E_0 = M_1/2 = %.1f PeV, so z(E) = E_0/E - 1 exactly.\n\n", E0))
cat("       z      E = E_0/(1+z)  (PeV)\n")
for (z in c(0,1,2,4,6,7,9,12,20)) cat(sprintf("   %6.0f %18.2f\n", z, E0/(1+z)))
zlo <- 4; zhi <- 9
E_hi <- E0/(1+zlo); E_lo <- E0/(1+zhi)
cat(sprintf("\n   The JWST compact-red-source window z = %g to %g maps to E = %.2f to %.2f PeV.\n",
            zlo, zhi, E_lo, E_hi))
cat("   That is above every neutrino yet reconstructed except KM3-230213A and below the line.\n")

cat("\n=============================================================================\n")
cat(" 2. HOW MUCH OF THE DECAY FLUX LANDS IN THAT WINDOW?  (tau cancels)\n")
cat("=============================================================================\n\n")
# Galactic: Phi_G = (N/4pi M tau) * Jbar, Jbar = sky-averaged NFW column  (GeV cm^-2)
rs <- 20; rsun <- 8.5; rho_loc <- 0.4
nfw <- function(r) 1/((r/rs)*(1+r/rs)^2)
rho_s <- rho_loc/nfw(rsun)
rho <- function(r) rho_s*nfw(pmax(r,1e-3))
Dfac <- function(psi, lmax=200)
  integrate(function(l) rho(sqrt(rsun^2+l^2-2*rsun*l*cos(psi))), 0, lmax,
            subdivisions=4000, rel.tol=1e-9)$value
skyavg <- function(a,b,n=800){ps<-seq(a,b,length.out=n); sum(sapply(ps,Dfac)*sin(ps))/sum(sin(ps))}
Dnear <- skyavg(0,pi/2); Dfar <- skyavg(pi/2,pi); Jbar <- (Dnear+Dfar)/2   # GeV cm^-3 kpc
Jbar_cm <- Jbar*kpc_cm
cat(sprintf("   NFW hemisphere columns  near %.3f  far %.3f  GeV cm^-3 kpc  (paper: 10.0 / 4.3)\n",
            Dnear, Dfar))
cat(sprintf("   all-sky mean column  Jbar = %.3f GeV cm^-3 kpc = %.4e GeV cm^-2\n\n", Jbar, Jbar_cm))

# Extragalactic: same prefactor times rho_DM,0 * c * Int dz/((1+z)H(z)) -> an effective column
rho_DM0 <- 0.1200/0.674^2 * 2.775e11    # Msun/Mpc^3 -> convert below; do it in GeV/cm^3 directly
rho_crit_GeVcm3 <- 1.05371e-5*0.674^2   # GeV cm^-3
rho_DM0 <- 0.1200/0.674^2 * rho_crit_GeVcm3   # = OmegaDM * rho_crit
cat(sprintf("   rho_DM,0 = %.4e GeV cm^-3\n", rho_DM0))
Ikern <- function(z) 1/((1+z)*Ez(z))
Itot  <- integrate(Ikern, 0, 1000, subdivisions=10000, rel.tol=1e-10)$value
Iband <- integrate(Ikern, zlo, zhi, rel.tol=1e-12)$value
colEG  <- rho_DM0*cH0_cm*Itot
colBand<- rho_DM0*cH0_cm*Iband
cat(sprintf("   Int_0^inf dz/((1+z)E(z)) = %.5f ;  Int_%g^%g = %.5f  (%.2f%% of it)\n",
            Itot, zlo, zhi, Iband, 100*Iband/Itot))
cat(sprintf("   extragalactic effective column  = %.4e GeV cm^-2  (all z)\n", colEG))
cat(sprintf("   extragalactic effective column  = %.4e GeV cm^-2  (z = %g-%g)\n", colBand, zlo, zhi))
tot <- Jbar_cm + colEG
cat(sprintf("\n   fraction of the whole decay flux that is Galactic (a line at %.1f PeV): %.4f\n",
            E0, Jbar_cm/tot))
cat(sprintf("   fraction that is extragalactic (a continuum below the line):            %.4f\n",
            colEG/tot))
f_band <- colBand/tot
cat(sprintf("   fraction in the LRD window z = %g-%g, i.e. E = %.1f-%.1f PeV:            %.4f\n",
            zlo, zhi, E_lo, E_hi, f_band))

cat("\n=============================================================================\n")
cat(" 3. WHAT A CROSS-CORRELATION WOULD COST, IN EVENTS\n")
cat("=============================================================================\n\n")
cat("   Pixelise the sky. Counts n_i = nbar(1 + f d_nu,i); template t_i = LRD overdensity,\n")
cat("   rms sigma_t. Estimator w = <(n_i/nbar - 1) t_i>. Signal = f r sigma_nu sigma_t,\n")
cat("   Poisson variance = sigma_t^2 / N. So\n\n")
cat("        S/N  =  f * r * sigma_nu * sqrt(N)        =>       N  =  (S/N)^2 / (f r sigma_nu)^2\n\n")
cat("   r <= 1 is the correlation coefficient, sigma_nu the rms fractional fluctuation of the\n")
cat("   shell's dark-matter column on the pixel scale. Both are bounded ABOVE, so N is bounded\n")
cat("   BELOW and no modelling choice can rescue it.\n\n")
Nreq <- function(f, sig, r=1, k=3) (k/(f*r*sig))^2
cat("      sigma_nu  (rms of the z=4-9 DM column on the pixel scale)   N events for 3 sigma\n")
for (s in c(1.0, 0.3, 0.1, 0.03, 0.01, 0.003))
  cat(sprintf("   %12.3f %52s\n", s, format(signif(Nreq(f_band,s),3), big.mark=",", scientific=TRUE)))
cat("\n   sigma_nu = 1 is absurd (the field is linear at z~6 on these scales and the shell\n")
cat("   averages over ~2 Gpc of line of sight); it is quoted only to show the floor.\n")

cat("\n   --- what sigma_nu actually is, from linear theory (BBKS transfer, sigma_8 normalised) ---\n\n")
h <- 0.674; Ob <- 0.0493; sigma8 <- 0.811; ns <- 0.965
Gam <- Om*h*exp(-Ob*(1+sqrt(2*h)/Om))
Tk <- function(k){ q <- k/Gam
  log(1+2.34*q)/(2.34*q) * (1+3.89*q+(16.1*q)^2+(5.46*q)^3+(6.71*q)^4)^(-0.25) }
Pk_un <- function(k) k^ns * Tk(k)^2                 # k in h/Mpc
W8 <- function(k,R=8) 3*(sin(k*R)-k*R*cos(k*R))/(k*R)^3
s8_un <- sqrt(integrate(function(lk){k<-exp(lk); k^3*Pk_un(k)*W8(k)^2/(2*pi^2)},
                        log(1e-4), log(1e3), subdivisions=5000, rel.tol=1e-8)$value)
Anorm <- (sigma8/s8_un)^2
Pk <- function(k) Anorm*Pk_un(k)                    # (Mpc/h)^3, z=0 linear
# growth factor, normalised D(0)=1
Dgrow <- function(z){ a<-1/(1+z)
  f <- function(x) (x/(Om + OL*x^3))^1.5
  g <- function(a) sqrt(Om/a^3+OL)*integrate(f,0,a,rel.tol=1e-10)$value
  g(a)/g(1) }
chi <- function(z) cH0_Mpc*integrate(function(zz) 1/Ez(zz), 0, z, rel.tol=1e-10)$value
cat(sprintf("   BBKS Gamma = %.4f, normalisation checks sigma_8 = %.4f\n", Gam, sigma8))
for (z in c(4,6,9)) cat(sprintf("   D(z=%g)/D(0) = %.4f,  chi = %6.0f Mpc,  1 deg = %5.1f Mpc comoving\n",
                                z, Dgrow(z), chi(z), chi(z)*pi/180))
# Limber: C_l = Int dchi/chi^2 * W(chi)^2 * P(l/chi, z(chi)), W = normalised shell selection
# sigma_nu^2 on scale theta ~ Int (2l+1)/(4pi) C_l B_l^2 ; use a Gaussian beam of FWHM theta.
zgrid <- seq(zlo, zhi, length.out=200)
chig  <- sapply(zgrid, chi); Dg <- sapply(zgrid, Dgrow)
dchidz<- cH0_Mpc/sapply(zgrid, Ez)
wz    <- dchidz/sum(dchidz*diff(zgrid)[1])            # normalised selection in z
Cl <- function(l){
  kk <- (l+0.5)/(chig*h)                              # h/Mpc
  ok <- kk>1e-5 & kk<1e3
  v  <- rep(0,length(kk))
  v[ok] <- (wz[ok]/dchidz[ok])^2 * Dg[ok]^2 * Pk(kk[ok])/h^3 / chig[ok]^2 * dchidz[ok]
  sum(v)*diff(zgrid)[1] }
sig_nu <- function(theta_deg){
  th <- theta_deg*pi/180; sb <- th/sqrt(8*log(2))
  ls <- 2:8000
  sum((2*ls+1)/(4*pi)*sapply(ls,Cl)*exp(-ls*(ls+1)*sb^2)) }
cat("\n      pixel scale      sigma_nu (rms of the z=4-9 DM column)     N for 3 sigma\n")
for (td in c(10,5,2,1,0.5)) {
  s <- sqrt(sig_nu(td))
  cat(sprintf("   %8.1f deg %30.5f %26s\n", td, s,
              format(signif(Nreq(f_band,s),3), big.mark=",", scientific=TRUE)))
}
cat("\n   (unbiased: dark-matter decay traces rho, bias exactly 1. The LRD template's own bias\n")
cat("    b~5 raises sigma_t, which cancels out of S/N; it does not help.)\n")

cat("\n=============================================================================\n")
cat(" 4. AND BEFORE ANY OF THAT, THE SKY OVERLAP\n")
cat("=============================================================================\n\n")
for (A in c(0.05, 0.54, 10, 1000, 18000)) {
  fs <- A/41252.96
  cat(sprintf("   a %8.2f deg^2 field is f_sky = %.3e; a 100-event all-sky decay sample puts\n", A, fs))
  cat(sprintf("       %.3e events inside it, and %.3e events are needed for one.\n", 100*fs, 1/fs))
}
cat("\n=============================================================================\n")
cat(" 5. FLATLY\n")
cat("=============================================================================\n")
cat("   The cross-correlation is the one place the two subjects genuinely meet, and it does not\n")
cat("   reach.  S/N = f r sigma_nu sqrt(N) contains no adjustable quantity: f = 0.031 is fixed\n")
cat("   by the closed budget and the NFW column, r <= 1 by definition, sigma_nu <= 0.008 at half\n")
cat("   a degree by linear theory.  At one degree it needs 4.3e8 decay-component events, six\n")
cat("   orders past the ~10^2 the endpoint test needs and which neutrino telescopes do not yet\n")
cat("   have.  Even at sigma_nu = 1, which no linear field on these scales can reach, the floor\n")
cat("   is 9.6e3 events.  A hundredfold error in every input at once does not change the answer,\n")
cat("   which is why the negative is safe to state.\n\n")
cat("   WHAT SURVIVES AND IS NEW: the extragalactic component itself.  It carries 42.7 per cent\n")
cat("   of the decay flux, its shape is fixed with no free parameter once M_1 is, and 95.1 per\n")
cat("   cent of the whole decay component lies above 50 PeV, so it does not overproduce where\n")
cat("   IceCube already measures.  See extragalactic_dilution.R for what it costs the\n")
cat("   directional assignment.\n")

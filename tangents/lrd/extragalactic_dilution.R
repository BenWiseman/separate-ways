# ==========================================================================================
# The framework's OWN decay model predicts an extragalactic component.  Section 3.2's
# directional-assignment numbers (57 events at 3 sigma, 159 at 5) assume "pure Galactic decay".
# The framework says that assumption is wrong by a factor of order two.  Price the correction.
# ==========================================================================================
M1 <- 491.6; E0 <- M1/2
H0 <- 67.4; Om <- 0.315; OL <- 0.685; c_kms <- 299792.458
Ez <- function(z) sqrt(Om*(1+z)^3+OL)
Mpc_cm <- 3.0856775814913673e24; kpc_cm <- Mpc_cm/1e3
cH0_cm <- c_kms/H0*Mpc_cm
rho_DM0 <- 0.1200/0.674^2 * 1.05371e-5*0.674^2

rs<-20; rsun<-8.5; rho_loc<-0.4
nfw<-function(r) 1/((r/rs)*(1+r/rs)^2); rho_s<-rho_loc/nfw(rsun)
rho<-function(r) rho_s*nfw(pmax(r,1e-3))
Dfac<-function(psi,lmax=200) integrate(function(l) rho(sqrt(rsun^2+l^2-2*rsun*l*cos(psi))),
        0,lmax,subdivisions=4000,rel.tol=1e-9)$value
skyavg<-function(a,b,n=800){ps<-seq(a,b,length.out=n); sum(sapply(ps,Dfac)*sin(ps))/sum(sin(ps))}
Dnear<-skyavg(0,pi/2); Dfar<-skyavg(pi/2,pi); Jbar<-(Dnear+Dfar)/2
p_gal <- Dnear/(Dnear+Dfar)
Jbar_cm <- Jbar*kpc_cm
cat(sprintf("  NFW hemisphere columns: near %.4f, far %.4f GeV cm^-3 kpc; p_near(pure Galactic) = %.4f\n",
            Dnear,Dfar,p_gal))
cat(sprintf("  all-sky mean Galactic column  = %.4e GeV cm^-2\n", Jbar_cm))

Ikern <- function(z) 1/((1+z)*Ez(z))
colEG <- function(z1,z2) rho_DM0*cH0_cm*integrate(Ikern,z1,z2,rel.tol=1e-12)$value
cat(sprintf("  all-z extragalactic column    = %.4e GeV cm^-2\n\n", colEG(0,1000)))

nev <- function(p,k) ceiling((k*0.5/(p-0.5))^2)
report <- function(tag, fEG) {
  p <- (1-fEG)*p_gal + fEG*0.5
  cat(sprintf("  %-52s f_EG = %.4f  p_near = %.4f   3sig %6d   5sig %6d\n",
              tag, fEG, p, nev(p,3), nev(p,5)))
}
cat("  CASE A: all decay-component events used, whatever their energy.\n")
report("    the paper's assumption (pure Galactic)", 0)
report("    what this framework actually predicts", colEG(0,1000)/(Jbar_cm+colEG(0,1000)))

cat("\n  CASE B: only events within one energy-resolution width of the line.\n")
cat("          A lognormal response of width s in ln E admits redshifts z < e^{k s} - 1.\n\n")
for (s in c(0.1,0.2,0.3)) for (k in c(1,2)) {
  zmax <- exp(k*s)-1
  fEG <- colEG(0,zmax)/(Jbar_cm+colEG(0,zmax))
  report(sprintf("    resolution %.0f%%, +/-%g sigma  (z < %.3f)", 100*s, k, zmax), fEG)
}
cat("\n  The dilution is unavoidable: the extragalactic component is the SAME line, from the\n")
cat("  same particle, at the same lifetime.  It cannot be switched off, only cut on energy,\n")
cat("  and cutting on energy throws away events the test is short of.\n")

cat("\n  --- the parameter-free spectral shape that comes with it ---\n\n")
cat("  dPhi/dE = K/(E H(z*)),  z* = E_0/E - 1,  E_0 = 245.8 PeV fixed by the abundance.\n")
cat("  Only the normalisation (1/tau) is free.  Normalise to the Galactic line integral:\n\n")
K <- rho_DM0*c_kms*1e5      # cgs-ish; ratios only, so the constant cancels
dPhidE <- function(E) { z<-E0/E-1; ifelse(z>0, 1/(E*Ez(z)), 0) }
gal_int <- Jbar_cm/(rho_DM0*cH0_cm)     # in the same units as Int dz/((1+z)E)
cat("        E (PeV)     z         E^2 dPhi/dE, relative to its value at E_0/2\n")
ref <- E0/2; refv <- ref^2*dPhidE(ref)
for (E in c(240,200,150,122.9,80,49.2,30,24.6,10,5,1))
  cat(sprintf("   %10.1f %8.2f %28.4f\n", E, E0/E-1, E^2*dPhidE(E)/refv))
cat("\n  The extragalactic continuum peaks in E^2 dPhi/dE near the line and falls as the\n")
cat("  expansion rate grows; it carries 42.7 per cent of the decay flux and its shape is a\n")
cat("  prediction with no free parameter beyond the overall lifetime.\n")

cat("\n  --- does the redshifted tail overproduce where IceCube already measures? ---\n")
# Flux per unit E is dPhi/dE ∝ 1/(E H(z)), z = E0/E - 1, so the tail runs to arbitrarily
# low energy. A referee will ask whether it floods the 10-100 TeV band IceCube has measured.
num <- integrate(function(E) dPhidE(E), 50, E0, rel.tol=1e-10)$value
den <- integrate(function(E) dPhidE(E), 1e-6, E0, rel.tol=1e-10, subdivisions=4000)$value
cat(sprintf("  fraction of the extragalactic decay component above 50 PeV : %.3f\n", num/den))
for (lo in c(0.01,0.1,1)) {
  f <- integrate(function(E) dPhidE(E), lo, E0, rel.tol=1e-10, subdivisions=4000)$value/den
  cat(sprintf("  fraction above %6.2f PeV : %.4f\n", lo, f))
}
cat("  The tail is steep because H(z) grows: most of the component stays near the line,\n")
cat("  so the prediction does not hide a flood at the energies already surveyed.\n")
fEG_all <- colEG(0,1000)/(Jbar_cm+colEG(0,1000))
cat(sprintf("\n  Combined with the Galactic line, which sits entirely at E_0:\n"))
cat(sprintf("  fraction of the WHOLE decay component above 50 PeV : %.3f\n",
            (1-fEG_all) + fEG_all*num/den))

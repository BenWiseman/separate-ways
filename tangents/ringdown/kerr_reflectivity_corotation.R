# The black-mirror horizon reflectivity is the GENERALISED Boltzmann factor,
#     |R| = exp( -|omega - m Omega_H| / (2 T_H) ),
# and the paper (and calc/tangents/blackmirror/boltzmann_reflectivity.py) dropped the
# co-rotation term m*Omega_H, using exp(-omega/2T_H). At zero spin Omega_H = 0 and the two
# agree, which is why the error hid. Away from zero they disagree in DIRECTION.
# Units M = 1, a = J/M^2, fundamental l = m = 2.
mw   <- function(a) 1.5251 - 1.1568*(1-a)^0.1292        # Berti-Cardoso-Will fit, as the paper uses
rp   <- function(a) 1 + sqrt(1-a^2)
TH   <- function(a) sqrt(1-a^2)/(4*pi*rp(a))
OmH  <- function(a) a/(2*rp(a))
m    <- 2
Rold <- function(a) exp(-mw(a)/(2*TH(a)))                       # what the paper used
Rnew <- function(a) exp(-abs(mw(a) - m*OmH(a))/(2*TH(a)))       # with co-rotation

cat("   a      M*omega    M*T_H     m*Omega_H   |R| paper    |R| correct   E returned\n")
for (a in c(0, 0.3, 0.7, 0.9, 0.95, 0.99)) {
  cat(sprintf(" %5.2f  %9.6f %9.6f %10.6f  %11.3e  %12.6f %12.6f\n",
      a, mw(a), TH(a), m*OmH(a), Rold(a), Rnew(a), Rnew(a)^2))
}
cat("\n  The paper says the reflectivity 'peaks at 9.8e-3 at zero spin and falls as spin rises'\n")
cat("  and that 'the horizon returns under 1e-4 of the incident energy at every spin'.\n")
cat("  With co-rotation it RISES with spin, and at a = 0.9 the horizon returns about 20 per\n")
cat("  cent of the incident energy.\n")

cat("\n  Does it still fail to reach the reported damping excess? |R| = 0.35 is what\n")
cat("  delta-tau_220 = 0.14 needs under the paper's own most generous damping model.\n\n")
need <- 0.35
cat("   a     |R| correct   reaches 0.35?   required T/T_H for |R| = 0.35\n")
for (a in c(0, 0.7, 0.9, 0.95)) {
  # |R| = exp(-|w-mOm|/(2T))  =>  T = |w-mOm| / (2 ln(1/|R|))
  Treq <- abs(mw(a) - m*OmH(a))/(2*log(1/need))
  cat(sprintf(" %5.2f  %12.6f   %-13s  %10.3f\n",
      a, Rnew(a), if (Rnew(a) >= need) "YES" else "no", Treq/TH(a)))
}
cat("\n  VERDICT: the paper's exclusion of the black mirror does not survive. At high spin the\n")
cat("  reflectivity fixed by T_H alone already exceeds what the damping excess asks for, and the\n")
cat("  required horizon is COOLER than Hawking rather than several times hotter. The claim that\n")
cat("  a T_H-fixed reflectivity cannot supply it is false. What remains true is narrower: a\n")
cat("  reflectivity magnitude is not a quasinormal-mode calculation, so this neither excludes\n")
cat("  the black mirror nor shows it explains the excess.\n")

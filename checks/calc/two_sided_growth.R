# Two-sided accretion multiplies the e-fold count, not the rate.
#
# A hole fed from both sheets takes Mdot_tot = (1+T) * Mdot_ours, where T is the fraction of the
# far sheet's infall the seam passes. Eddington-limited growth is exponential with e-folding time
# t_S/(1+T), so over a FIXED time window the number of e-folds scales linearly in 1+T. The fold's
# own time reversal forces T = 0 unless matter crosses inside the contact radius, which is what
# checks/calc/two_sheet_accretion.py settles. That is the whole content: a factor
# two in the rate is a factor two in the exponent, which is why it matters at all.
#
# Everything below is computed. The planted failure is at the end.

c_cgs   <- 2.99792458e10; G <- 6.67430e-8
sigT    <- 6.6524587e-25; mp <- 1.67262192e-24
Myr     <- 3.1557e13          # s
Gyr     <- 1000*Myr

salpeter <- function(eps) (eps/(1-eps)) * sigT*c_cgs/(4*pi*G*mp)

# --- flat LCDM age, Planck 2018 ---
H0 <- 67.36; Om <- 0.3153; Ol <- 1-Om
H0_s <- H0*1e5/3.0857e24                      # s^-1
age <- function(z) integrate(function(zz) 1/((1+zz)*H0_s*sqrt(Om*(1+zz)^3+Ol)),
                             z, Inf, rel.tol=1e-10)$value

z_hi <- 20; z_lo <- 7
dt <- age(z_lo) - age(z_hi)

cat(sprintf("\n  age(z=%d) = %.3f Gyr,  age(z=%d) = %.3f Gyr,  window = %.4f Gyr\n",
            z_hi, age(z_hi)/Gyr, z_lo, age(z_lo)/Gyr, dt/Gyr))

for (eps in c(0.057, 0.1, 0.2)) {
  tS <- salpeter(eps)
  n1 <- dt/tS
  cat(sprintf("  eps = %.3f : t_Salpeter = %5.1f Myr, one-sided e-folds = %6.2f, two-sided (T=1) = %6.2f, mass ratio = %.3g\n",
              eps, tS/Myr, n1, 2*n1, exp(n1)))
}

tS <- salpeter(0.1); n_efold <- dt/tS
cat(sprintf("\n  HEADLINE at eps = 0.1: %.2f e-folds from z=20 to z=7 with the seam closed.\n", n_efold))
cat(sprintf("  A seed of 100 Msun reaches %.3g Msun one-sided and %.3g Msun two-sided.\n",
            100*exp(n_efold), 100*exp(2*n_efold)))
cat(sprintf("  Equivalently: two-sided growth reaches the one-sided endpoint in half the time,\n"))
cat(sprintf("  i.e. by z = %.1f rather than z = %d.\n",
            {f <- function(z) (age(z)-age(z_hi)) - dt/2; uniroot(f, c(7,20))$root}, z_lo))

# --- how big must T be to matter? ---
cat("\n  transmitted fraction T needed to gain a given mass factor over the same window:\n")
for (gain in c(10, 100, 1000)) {
  k <- 1 + log(gain)/n_efold
  cat(sprintf("     x%-5g  needs T = %.3f  (i.e. the seam passing %.1f%% of the far sheet's infall)\n",
              gain, k, 100*(k-1)))
}

# --- PLANTED FAILURE: if the multiplier hit the RATE-as-prefactor rather than the exponent,
#     M would go as (1+T)exp(t/tS) and the gain would be 1+T itself, not exp(T n).
wrong <- 2
right <- exp((2-1)*n_efold)
cat(sprintf("\n  CHECK (planted): prefactor reading gives a x%.0f gain, exponent reading gives x%.3g.\n",
            wrong, right))
stopifnot(right/wrong > 1e4)
cat("  Ratio is 1e4 or more, so the two readings are not confusable and the exponent one is used.\n")

cat("\n=== flatly ===\n")
cat(sprintf("
  Over the window the little red dots occupy, z = 20 to z = 7, an Eddington-limited
  hole gets %.1f e-folds at radiative efficiency 0.1. Feeding it from both sheets
  does not add %.0f per cent to its mass, it doubles the exponent: the same window
  buys %.1f e-folds, and a %.3g times larger hole. That is the only place in this
  paper where the second sheet changes an astrophysical number by a large factor,
  and it is also the softest, because T is not derived here. Fold invariance fixes the
  other sheet's accretion rate to equal ours, so T is a transmitted fraction and not a
  rate ratio, but what sets it is the contact geometry and that is not computed.

  What does not depend on T's value is the shape of the statement. Any T > 0
  at all raises the exponent, so the two-sided class is systematically overmassive
  at fixed seed and fixed epoch, and a mass gain of a hundred needs only T = %.2f.
  The claim to make is that direction, not the factor.\n",
  n_efold, 100*(2-1), 2*n_efold, exp(n_efold), log(100)/n_efold))

# ---------------------------------------------------------------------------
# WHICH SALPETER TIME. This is a real fork and it moved the headline by 4x.
# An earlier note carried 12.93 e-folds and a 4e5 mass factor. That comes from
# t_S = 45 Myr at eps = 0.1, i.e. dropping the (1-eps) in the denominator.
#
# For BLACK HOLE MASS growth the (1-eps) belongs there: the hole gains
# Mdot_BH = (1-eps) Mdot_inflow while the luminosity L = eps Mdot_inflow c^2
# is what the Eddington limit caps. So
#     t_S = M/Mdot_BH = (eps/(1-eps)) * sigma_T c / (4 pi G m_p),
# which is 50.1 Myr at eps = 0.1, not 45. The 45 Myr form is the common
# shorthand and it overstates the growth.
cat("\n=== which Salpeter time ===\n")
t_short <- 0.1 * sigT*c_cgs/(4*pi*G*mp)       # the (1-eps)-free shorthand
t_full  <- salpeter(0.1)
cat(sprintf("  shorthand t_S = %.1f Myr -> %.2f e-folds -> mass factor %.3g\n",
            t_short/Myr, dt/t_short, exp(dt/t_short)))
cat(sprintf("  correct   t_S = %.1f Myr -> %.2f e-folds -> mass factor %.3g\n",
            t_full/Myr,  dt/t_full,  exp(dt/t_full)))
stopifnot(abs(dt/t_short - 12.9) < 0.1)        # reproduces the note's number
stopifnot(abs(dt/t_full  - 11.6) < 0.1)        # and the corrected one
cat("  Both reproduced, so the discrepancy is the (1-eps) and nothing else.\n")

# ---------------------------------------------------------------------------
# 5. Turning T round: the observations bound it, rather than it being chosen.
#
# Everything above runs forwards, from a chosen T to a mass gain. The gain is also an
# OBSERVABLE, because the far sheet feeds the hole and not our galaxy: two-sided growth raises
# M_BH at fixed stellar mass, so the mass gain is the overmassiveness relative to a one-sided
# hole with the same seed and the same host.
#
# JWST's overmassive early holes sit roughly one to two orders above the local relation. Ask what
# T that is, and what a seam fully open to matter would have predicted.

efolds1 <- (tz7 <- age(7)) - (tz20 <- age(20))
efolds1 <- efolds1 / tS                      # one-sided e-folds over the window
k_for   <- function(x) 1 + log(x)/efolds1
gain_at <- function(k) exp((k - 1) * efolds1)

cat("\n=== 5. what the observed overmassiveness asks of T ===\n\n")
cat("   overmassive by   implied T   seam passing\n")
for (x in c(3, 10, 30, 100, 300)) {
  k <- k_for(x)
  cat(sprintf("   %10.0f x  %10.3f   %5.1f per cent of the far sheet's infall\n",
              x, k - 1, 100*(k-1)))
}
cat(sprintf("\n   and a seam fully open to matter, T = 1, predicts a gain of %.1e,\n", gain_at(2)))
cat(sprintf("   which overshoots a hundredfold overmassiveness by a factor of %.0e.\n",
            gain_at(2)/100))
stopifnot(abs(k_for(10) - 1.198) < 2e-3, abs(k_for(100) - 1.396) < 2e-3, gain_at(2) > 1e4)

# --- 6. the shape of the population, which is what distinguishes this from the menu ---
#
# Heavy seeds, super-Eddington accretion, primordial seeds and mass-measurement bias all move
# the WHOLE population. This does not: only the two-sided class is affected, and T is a property
# of the seam rather than of the individual hole, so at Lambda M^2 ~ 0 every two-sided hole gets
# the SAME multiplier. The offset in dex is T times the one-sided e-fold count over log(10), so
# the prediction is two loci at a fixed separation rather than one locus with extra scatter.

dex <- function(T) T * efolds1 / log(10)
cat("\n=== 6. offset of the two-sided locus, at fixed seed and fixed epoch ===\n\n")
cat(sprintf("   one-sided e-folds over the window: %.2f\n", efolds1))
for (T in c(0.1, 0.198, 0.396, 1.0))
  cat(sprintf("   T = %.3f  ->  %.2f dex above the one-sided locus\n", T, dex(T)))
stopifnot(abs(dex(0.198) - 1.00) < 0.02, abs(dex(0.396) - 2.00) < 0.02)
cat("\n   So the two classes are separated, not blended: the same T displaces every two-sided\n")
cat("   hole of the same growth history by the same factor. A mechanism that shifts the mean\n")
cat("   cannot make two loci, and a complete sample at fixed stellar mass tells them apart.\n")


cat("\n   So if these objects were two-sided, the data would fix T near 0.2 to 0.4 and\n")
cat("   EXCLUDE the fully open seam by three to four orders of magnitude. The free parameter\n")
cat("   of 6.1 is then bounded by observation rather than chosen, which is the opposite of how\n")
cat("   it enters above. What the argument assumes is that the far sheet feeds the hole and not\n")
cat("   our galaxy, so the host's stellar mass is unaffected and the mass gain IS the\n")
cat("   overmassiveness; that is what a second exterior means.\n")
cat("\n   It stays conditional on there being sub-maximal two-sided holes at all, which 6.2\n")
cat("   leaves open outside the Schwarzschild-de Sitter family and does not supply.\n")

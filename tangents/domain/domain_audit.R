# A.14 says the appendix belongs to the maximally extended, equilibrium case. That may be
# too blunt: the claims do not all lean on the same thing. Audit them by DEPENDENCE, and
# test the ones that matter rather than reasoning about them.
#
# The decisive test: the reciprocal law's proof uses only P_perp^2 = 1 and positivity. If
# that is right it cannot care whether the state is thermal, so it should survive an
# arbitrarily non-KMS state. Build one and check.

set.seed(41); LMAX <- 200
Pl <- function(l,x){ if(l==0) return(1); if(l==1) return(x)
  p0<-1; p1<-x; for(k in 1:(l-1)){p2<-((2*k+1)*x*p1-k*p0)/(k+1); p0<-p1; p1<-p2}; p1 }
Wsum <- function(Rl, g, sgn) { l <- 0:LMAX
  sum((2*l+1)/(4*pi) * Rl * sgn^l * sapply(l, function(k) Pl(k, cos(g)))) }

cat("=== 1. does the reciprocal law need thermality? Build states that are not thermal.\n\n")
beta <- 1
states <- list(
  "thermal, KMS"                = 1/(exp(beta*(0:LMAX)/8)-1+1e-12),
  "Unruh-like, flux-carrying"   = (1/(exp(beta*(0:LMAX)/8)-1+1e-12))*(1+0.4*(0:LMAX)/(1+(0:LMAX))),
  "wildly non-thermal"          = abs(rnorm(LMAX+1, 1, 1))*exp(-(0:LMAX)/25),
  "non-monotone, oscillating"   = (1.5+cos((0:LMAX)/3))*exp(-(0:LMAX)/30),
  "power law, no temperature"   = (1+(0:LMAX))^(-1.7)
)
cat("        state                          max |R(x,y) R(x,Py) - 1| over 5 angles\n")
for (nm in names(states)) {
  Rl <- states[[nm]]
  e <- max(sapply(c(25,60,95,130,165)*pi/180, function(g) {
      Rxy <- Wsum(Rl,g,-1)/Wsum(Rl,g,+1)          # P_perp inserts (-1)^l
      Rxp <- Wsum(Rl,pi-g,-1)/Wsum(Rl,pi-g,+1)    # the image point
      abs(Rxy*Rxp - 1) }))
  cat(sprintf("   %-32s %.3e\n", nm, e))
}
cat("\n  Machine precision for every state, thermal or not. The reciprocal law does not\n")
cat("  care about temperature, flux, or whether a temperature exists at all. So it is\n")
cat("  NOT restricted to the equilibrium case and A.14's blanket caveat overstates it.\n")

cat("\n=== 2. what each claim actually leans on\n\n")
cat("        claim                              needs KMS   needs B   needs two sheets\n")
rows <- list(
 c("reciprocal law R R' = 1",              "no",  "no",  "no"),
 c("ninety-degree equality",               "no",  "no",  "no"),
 c("log R has odd multipoles only",        "no",  "no",  "no"),
 c("reduction to W_B",                     "no",  "yes", "no"),
 c("alpha^2 = 1 fixes the temperature",    "YES", "no",  "yes"),
 c("Hadamard forces a transparent seam",   "no",  "no",  "yes"),
 c("mass-independent reflectivity",        "no",  "no",  "no"),
 c("F is the unique meeting region",       "no",  "no",  "yes"),
 c("singularity consistency condition",    "no",  "no",  "yes"))
for (r in rows) cat(sprintf("   %-36s %9s %9s %16s\n", r[1], r[2], r[3], r[4]))
cat("\n  Only ONE claim needs KMS. FOUR need the two-sided structure, the KMS one among\n")
cat("  them -- an earlier version of this line said three and dropped it, and the table\n")
cat("  above is what is right. One needs B.\n")
cat("  So evaporation, which breaks KMS, threatens exactly one line of the appendix.\n")

cat("\n=== 3. and how badly is KMS broken for a real black hole?\n\n")
G <- 6.67430e-11; c_ <- 2.99792458e8; hbar <- 1.054571817e-34; Msun <- 1.98847e30
lP <- sqrt(hbar*G/c_^3); tP <- lP/c_; mP <- sqrt(hbar*c_/G)
cat("   The departure from equilibrium is of order 1/(kappa t_evap). With kappa ~ 1/(4M)\n")
cat("   and t_evap ~ M^3 in Planck units, that is ~ M^-2, M in Planck masses.\n\n")
cat("        M (Msun)      M / m_Planck        1/(kappa t_evap) ~ M^-2\n")
for (M in c(3, 10, 65, 1e6, 1e9)) {
  Mp <- M*Msun/mP
  cat(sprintf("   %12.3g %18.3e %28.3e\n", M, Mp, 1/Mp^2))
}
cat("\n  Between 1e-77 and 1e-93. So for the ONE claim that needs KMS, the departure is\n")
cat("  smaller than any number elsewhere in the paper by fifty orders of magnitude.\n")

cat("\n=== 4. but small is not the same as smooth, and last night's lesson applies\n\n")
cat("  The singularity consistency condition broke DISCONTINUOUSLY: any epsilon gave an\n")
cat("  empty solution space, not a small one. So the question for the temperature claim\n")
cat("  is which behaviour it has. Test it: perturb the state away from thermal and\n")
cat("  measure how far W(t - i beta) - W(t) moves.\n\n")
H <- 1; bb <- 2*pi/H
cat("  A first attempt used sech^2(Ht/2) as the admixture. That function is ITSELF\n")
cat("  beta-periodic (cosh(x - i pi) = -cosh x), so it broke nothing and the test was\n")
cat("  vacuous: the residual sat at 4e-15 for every admixture including 0.5. Replaced\n")
cat("  with a term carrying no imaginary period at all.\n\n")
W <- function(t, d) 1/sinh(H*t/2)^2 + d*exp(-t^2)        # exp(-t^2) has no imaginary period
cat("        d           max |W(t - i beta) - W(t)|      ratio to previous d\n")
prev <- NA
for (d in c(0, 1e-8, 1e-4, 1e-2, 0.5)) {
  v <- max(sapply(c(0.8,1.5,2.6,3.4), function(t) Mod(W(t-1i*bb,d) - W(t,d))))
  cat(sprintf("   %8.0e %28.3e %22s\n", d, v, if (is.na(prev)||prev==0) "-" else sprintf("%.1f", v/prev)))
  prev <- v
}
cat("\n  Exactly linear in the admixture: the ratios track d. So the behaviour is\n")
cat("  CONTINUOUS through zero, unlike the singularity condition which collapsed to an\n")
cat("  empty solution space at any epsilon. That qualitative difference is the result.\n")
cat("\n  The COEFFICIENT is not established and the numbers above overstate it. A Gaussian\n")
cat("  blows up under imaginary time translation: exp(-(t-i beta)^2) carries exp(beta^2),\n")
cat(sprintf("  which is %.1e here, so the large absolute values are the test function and not\n", exp(bb^2)))
cat("  physics. Getting the physical coefficient needs the actual non-equilibrium\n")
cat("  correction to the state, which is not computed here.\n")
cat("\n  What survives the ignorance: even carrying that inflated coefficient, a departure\n")
cat(sprintf("  of 1e-77 gives a KMS violation of order %.0e. The conclusion is robust across\n", 1e-77*exp(bb^2)))
cat("  more than fifty orders of magnitude of coefficient, so the controlled-approximation\n")
cat("  claim does not rest on knowing it.\n")

# Is A.18's n = 1/2 distinctive, or would decoherence at the bang produce it anyway?
#
# The obvious referee question. A.18 derives n = 1/2 at the contact from Theta-invariance.
# Dephasing during a Landau-Zener sweep also moves populations around, and the quantum Zeno
# effect freezes them. If a strongly dephased sweep gave 1/2 on its own, the fold's signature
# would not be distinctive. Compute it, in BOTH dephasing bases, since the answer differs.
#
# H = (1/2)(eps(t) sz + Delta sx), eps = gamma t, started in the lower state at t = -T.
run <- function(gamma, Delta, Gam, T=60, n=240000, basis="z") {
  h <- 2*T/n; t <- -T; r <- c(0,0,1)
  f <- function(t,r) {
    eps <- gamma*t
    if (basis=="z")                       # dephasing in the DIABATIC basis
      c(-eps*r[2]-Gam*r[1], eps*r[1]-Delta*r[3]-Gam*r[2], Delta*r[2])
    else {                                # dephasing along the instantaneous field (ADIABATIC)
      B <- c(Delta,0,eps); Bn <- B/sqrt(sum(B^2)); perp <- r - sum(r*Bn)*Bn
      c(-eps*r[2], eps*r[1]-Delta*r[3], Delta*r[2]) - Gam*perp
    }
  }
  for (i in seq_len(n)) {
    k1<-f(t,r); k2<-f(t+h/2,r+h/2*k1); k3<-f(t+h/2,r+h/2*k2); k4<-f(t+h,r+h*k3)
    r <- r + h/6*(k1+2*k2+2*k3+k4); t <- t+h
  }
  (1-r[3])/2
}
cat("=== 1. coherent limit reproduces Landau-Zener ===\n")
cat("  the integrator returns the ADIABATIC probability, so compare with 1 - exp(-pi D^2/2g):\n")
cat("     Delta   numeric (Gam=0)   1 - exp(-pi D^2/2g)\n")
for (D in c(0.4,0.8,1.2))
  cat(sprintf("   %7.2f %17.6f %21.6f\n", D, run(1,D,0), 1-exp(-pi*D^2/2)))
cat("  (residual is finite-T truncation of the sweep, not method error.)\n")

cat("\n=== 2. dephasing in the DIABATIC basis: Zeno freezes the transition, n -> 0 ===\n")
cat("     Delta     Gam=0    Gam=10    Gam=50   Gam=200\n")
for (D in c(0.4,0.8,1.2)) { cat(sprintf("   %7.2f",D))
  for (G in c(0,10,50,200)) cat(sprintf(" %9.5f", run(1,D,G))); cat("\n") }

cat("\n=== 3. dephasing in the ADIABATIC basis: the state follows, n -> 1 ===\n")
cat("     Delta     Gam=0     Gam=1     Gam=10    Gam=50\n")
for (D in c(0.4,0.8,1.2)) { cat(sprintf("   %7.2f",D))
  for (G in c(0,1,10,50)) cat(sprintf(" %9.5f", run(1,D,G,basis="a"))); cat("\n") }

cat("\n  VERDICT, and it runs the paper's way. The two strong-decoherence limits go to OPPOSITE\n")
cat("  ends: diabatic-basis dephasing freezes the populations and drives n to 0, adiabatic-basis\n")
cat("  dephasing enforces following and drives n to 1. One half is the one value NEITHER limit\n")
cat("  reaches, and it is not approached from either side as Gam grows. So decoherence does not\n")
cat("  mimic A.18's result, and a referee asking whether n = 1/2 is just dephasing has an answer.\n")
cat("\n  What this does NOT show. It is the same two-level reduction A.18 works in, with a\n")
cat("  phenomenological Lindblad dephasing rate rather than one derived from the cosmological\n")
cat("  environment, and with the sweep truncated at finite time. It rules out the easy\n")
cat("  confound; it does not prove no environment can produce one half by some other route.\n")

cat("\n  ==========================================================================\n")
cat("  CORRECTION 2026-09-21, third pre-submission pass. THE CONCLUSION ABOVE IS WRONG.\n")
cat("  ==========================================================================\n")
cat("  Every number above uses T = 60 (see the default in run()). That window is far too\n")
cat("  short at strong dephasing, and the result is not converged in it. Re-running the same\n")
cat("  Bloch equation at longer windows, with the step bounded by the largest rate:\n\n")
bl <- function(Delta, Gam, v, T, pps=0.1) {
  nstep <- ceiling(2*T*max(2*Delta, v*T, Gam)/pps)
  dt <- 2*T/nstep; t <- -T; r <- c(0,0,-1)
  for (k in 1:nstep) {
    f <- function(tt, rr) c(-v*tt*rr[2] - Gam*rr[1],
                             v*tt*rr[1] - 2*Delta*rr[3] - Gam*rr[2], 2*Delta*rr[2])
    k1<-f(t,r); k2<-f(t+dt/2,r+dt/2*k1); k3<-f(t+dt/2,r+dt/2*k2); k4<-f(t+dt,r+dt*k3)
    r <- r + dt/6*(k1+2*k2+2*k3+k4); t <- t+dt }
  (1 + r[3])/2
}
cat("      Delta   Gamma      T=60       T=300      T=900\n")
for (D in c(0.4, 2.0))
  cat(sprintf("   %8.2f %7.0f %10.6f %10.6f %10.6f\n", D, 200,
              bl(D,200,1,60), bl(D,200,1,300), bl(D,200,1,900)))
cat("\n  At Delta = 0.4 the answer climbs from 0.156 to 0.411 and is still rising: the 0.045\n")
cat("  quoted above is a T = 60 artefact. At Delta = 2 it is 0.5000000, EXACTLY the value the\n")
cat("  text claims neither limit reaches.\n")
cat("\n  So strong diabatic-basis dephasing DOES reach one half, and the easy confound is NOT\n")
cat("  excluded. This agrees with the known strong-dephasing Landau-Zener result (Saito and\n")
cat("  Kayanuma, cond-mat/0111420, eq. 25), which is half the coherent transition probability\n")
cat("  and tends to one half as the coupling grows.\n")
cat("\n  What survives: the two basis limits still DIFFER, diabatic dephasing giving one half\n")
cat("  and adiabatic dephasing driving toward one. What does not: the claim that one half is\n")
cat("  unreachable by an environment. B.3 is rewritten to say so.\n")

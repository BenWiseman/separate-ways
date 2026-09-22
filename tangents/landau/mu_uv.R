# mu* moves with p, so what Theta leaves is a phase FUNCTION mu(p), not one number. The
# question that matters: does the ultraviolet condition 5.1 already uses pick that function?
# Fixed mu should leave the half-sweep p^-4 boundary tail. Only the minimising mu*(p) can
# cancel it down to the Gaussian. If so, Hadamard forces mu(p) and the import closes.
#
# Precision note stated up front: at p the Gaussian is exp(-pi p^2) while the power law is
# ~p^-4/16, so the cancellation needed is ~1e-3 at p=2, 2e-6 at p=2.5, 6e-10 at p=3 and
# 1e-18 at p=4. Double precision and RK4 cannot resolve beyond about p=3, and the table
# below stops where the integrator stops being trustworthy rather than where the claim stops.

gamma <- 1
prop <- function(p, eta0, eta1, nstep) {
  h <- (eta1-eta0)/nstep
  f <- function(eta,M) rbind(-1i*(gamma*eta*M[1,]+p*M[2,]), -1i*(p*M[1,]-gamma*eta*M[2,]))
  M <- diag(2)+0i; eta <- eta0
  for (s in seq_len(nstep)) {
    k1<-f(eta,M);k2<-f(eta+h/2,M+h/2*k1);k3<-f(eta+h/2,M+h/2*k2);k4<-f(eta+h,M+h*k3)
    M<-M+h/6*(k1+2*k2+2*k3+k4); eta<-eta+h }
  M
}
vecU <- function(eta,p){ E<-sqrt((gamma*eta)^2+p^2); U<-c(p+0i,-gamma*eta+E); U/sqrt(sum(Mod(U)^2)) }

Tf <- 14; N <- 900000
ps <- c(1.0,1.25,1.5,1.75,2.0,2.25,2.5,2.75,3.0)
nmin <- nfix <- npi <- numeric(length(ps))
for (i in seq_along(ps)) {
  p <- ps[i]; Uh <- prop(p,0,Tf,N); uT <- vecU(Tf,p)
  a1 <- sum(Conj(uT)*(Uh%*%c(1,0))); a2 <- sum(Conj(uT)*(Uh%*%c(0,1)))
  nmin[i] <- (Mod(a1)-Mod(a2))^2/2          # minimum over mu, attained at mu*(p)
  nfix[i] <- Mod(a1+exp(1i*2.0)*a2)^2/2     # one FIXED phase, mu = 2.0 rad
  npi[i]  <- Mod(a1-a2)^2/2                 # mu = pi, the bang-adiabatic state
}
cat("     p    min over mu      fixed mu=2.0      mu=pi (bang-ad)   exp(-pi p^2)\n")
for (i in seq_along(ps))
  cat(sprintf("  %5.2f %14.4e %16.4e %18.4e %14.4e\n",
      ps[i], nmin[i], nfix[i], npi[i], exp(-pi*ps[i]^2)))

sl <- function(v) diff(log(v))/diff(log(ps))
cat("\n  log-slopes d ln n / d ln p (p^-4 predicts -4; Gaussian predicts -2 pi p^2):\n")
cat(sprintf("    min over mu : %s\n", paste(sprintf("%7.2f", sl(nmin)), collapse="")))
cat(sprintf("    fixed mu=2.0: %s\n", paste(sprintf("%7.2f", sl(nfix)), collapse="")))
cat(sprintf("    mu = pi     : %s\n", paste(sprintf("%7.2f", sl(npi)),  collapse="")))
cat(sprintf("    Gaussian    : %s\n", paste(sprintf("%7.2f", -2*pi*ps[-1]^2), collapse="")))
cat("\n  ratio of min-over-mu to exp(-pi p^2), which should be 1/4 if the cancellation is exact:\n")
cat(sprintf("    %s\n", paste(sprintf("%8.4f", nmin/exp(-pi*ps^2)), collapse="")))

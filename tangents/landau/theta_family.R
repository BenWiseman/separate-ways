# The fold on this mode problem. H(eta) = [[g eta, p],[p, -g eta]] is real and symmetric, and
# H(-eta) = sigma_x H(eta) sigma_x, so the antilinear map  Theta: psi(eta) -> sigma_x psi*(-eta)
# is a symmetry with Theta^2 = 1. That is the fold, written for one mode.
#
# A Theta-invariant state obeys psi(0) = e^{i d} sigma_x psi*(0), which forces
# |psi_1(0)| = |psi_2(0)| = 1/sqrt2: the mode is exactly half in each component AT THE BANG.
# The ray is then psi(0) ~ (1, e^{i mu})/sqrt2, one real parameter. Map the family.
# Evolution is linear, so integrate the two basis columns once and every mu is free.

gamma <- 1
prop <- function(p, eta0, eta1, nstep) {      # 2x2 propagator, columns = images of e1,e2
  h <- (eta1-eta0)/nstep
  f <- function(eta,M) rbind(-1i*(gamma*eta*M[1,] + p*M[2,]),
                             -1i*(p*M[1,] - gamma*eta*M[2,]))
  M <- diag(2)+0i; eta <- eta0
  for (s in seq_len(nstep)) {
    k1<-f(eta,M);k2<-f(eta+h/2,M+h/2*k1);k3<-f(eta+h/2,M+h/2*k2);k4<-f(eta+h,M+h*k3)
    M<-M+h/6*(k1+2*k2+2*k3+k4); eta<-eta+h }
  M
}
vecs <- function(eta,p){ E<-sqrt((gamma*eta)^2+p^2)
  L<-c(p+0i,-gamma*eta-E); U<-c(p+0i,-gamma*eta+E)
  list(L=L/sqrt(sum(Mod(L)^2)),U=U/sqrt(sum(Mod(U)^2))) }

Tf<-16; N<-200000
mus <- seq(0,2*pi,length.out=2001)
cat("    p     P=|beta|^2   half-angle n    min over mu    max over mu   mu* at half-angle   n(mu=pi)\n")
for (p in c(0.3,0.5,0.7,0.9,1.1)) {
  vin<-vecs(-Tf,p); vout<-vecs(Tf,p)
  Ufull <- prop(p,-Tf,Tf,N)
  P  <- Mod(sum(Conj(vout$U)*(Ufull%*%vin$L)))^2
  nh <- (1-sqrt(1-P))/2
  Uh <- prop(p,0,Tf,N/2)
  c1 <- Uh%*%c(1,0); c2 <- Uh%*%c(0,1)
  a1 <- sum(Conj(vout$U)*c1); a2 <- sum(Conj(vout$U)*c2)
  ns <- Mod(a1 + exp(1i*mus)*a2)^2/2
  istar <- which.min(abs(ns-nh))
  nPi <- Mod(a1 - a2)^2/2
  cat(sprintf("  %4.2f %12.7f %14.7f %14.7f %14.7f %14.4f %14.3e\n",
      p, P, nh, min(ns), max(ns), mus[istar], nPi))
}
cat("\n  Reading: Theta fixes the magnitudes at the bang and leaves the relative phase mu free.\n")
cat("  Whether the half-angle occupation is the MINIMUM of that family decides whether BFT's\n")
cat("  'minimise the occupation' is a minimum over a family the FOLD defines, rather than a\n")
cat("  prescription imposed on all states.\n")

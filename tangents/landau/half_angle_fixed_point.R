# A.1 says fold-invariance is "the exchange of in- and out-regions". A.5 calls BFT's state
# the "half-angle vacuum". If those two facts are the same fact, the half angle is the FIXED
# POINT of the exchange and 5.1's prescription is forced, not imported. Test it against the
# real evolution operator, with the Stokes phase included rather than assumed away.

gamma <- 1
rk4v <- function(p, eta0, eta1, psi0, nstep) {
  h <- (eta1-eta0)/nstep
  f <- function(eta, psi) c(-1i*(gamma*eta*psi[1] + p*psi[2]),
                            -1i*(p*psi[1] - gamma*eta*psi[2]))
  psi <- psi0; eta <- eta0
  for (s in seq_len(nstep)) {
    k1<-f(eta,psi); k2<-f(eta+h/2,psi+h/2*k1); k3<-f(eta+h/2,psi+h/2*k2); k4<-f(eta+h,psi+h*k3)
    psi <- psi + h/6*(k1+2*k2+2*k3+k4); eta <- eta+h
  }
  psi
}
vecs <- function(eta,p){ E<-sqrt((gamma*eta)^2+p^2)
  L<-c(p+0i,-gamma*eta-E); U<-c(p+0i,-gamma*eta+E)
  list(L=L/sqrt(sum(Mod(L)^2)), U=U/sqrt(sum(Mod(U)^2))) }

Tf <- 16; N <- 240000
cat("   p      P=|beta|^2    half-angle n      equal-occupation n     |diff|     phases tried\n")
for (p in c(0.3,0.5,0.7,0.9,1.1)) {
  vin <- vecs(-Tf,p); vout <- vecs(Tf,p)
  Lin_T <- rk4v(p,-Tf,Tf,vin$L,N)      # in-vacuum pushed to +T
  Uin_T <- rk4v(p,-Tf,Tf,vin$U,N)      # in-excited pushed to +T
  P  <- Mod(sum(Conj(vout$U)*Lin_T))^2
  nh <- (1-sqrt(1-P))/2                # BFT half-angle prediction

  # general vacuum at +T:  |v> = cos th |L_out> + e^{i ph} sin th |U_out>
  # n_out = sin^2 th ;  n_in = |<Uin_T | v>|^2 .  Solve n_in = n_out over (th, ph).
  best <- NULL
  for (ph in seq(0, 2*pi, length.out = 721)) {
    g <- function(th) { v <- cos(th)*vout$L + exp(1i*ph)*sin(th)*vout$U
                        Mod(sum(Conj(Uin_T)*v))^2 - sin(th)^2 }
    lo <- g(1e-6); hi <- g(pi/2-1e-6)
    if (is.finite(lo) && is.finite(hi) && lo*hi < 0) {
      th <- uniroot(g, c(1e-6, pi/2-1e-6), tol=1e-12)$root
      n  <- sin(th)^2
      if (is.null(best) || abs(n-nh) < abs(best$n-nh)) best <- list(n=n, ph=ph)
    }
  }
  # how much does the solution move as the phase is scanned?
  ns <- c()
  for (ph in seq(0, 2*pi, length.out = 181)) {
    g <- function(th) { v <- cos(th)*vout$L + exp(1i*ph)*sin(th)*vout$U
                        Mod(sum(Conj(Uin_T)*v))^2 - sin(th)^2 }
    lo<-g(1e-6); hi<-g(pi/2-1e-6)
    if (is.finite(lo)&&is.finite(hi)&&lo*hi<0) ns <- c(ns, sin(uniroot(g,c(1e-6,pi/2-1e-6),tol=1e-12)$root)^2)
  }
  cat(sprintf("  %4.2f %12.7f %16.7f %20.7f %11.2e   spread over phase: %.2e\n",
      p, P, nh, if(is.null(best)) NA else best$n,
      if(is.null(best)) NA else abs(best$n-nh), if(length(ns)) diff(range(ns)) else NA))
}
cat("\n  Reading: if 'equal-occupation n' equals the half-angle n and the spread over the\n")
cat("  Stokes phase is zero, then the state that looks equally excited from the in- and\n")
cat("  out-regions IS BFT's state, uniquely, and 5.1's prescription is forced by the fold.\n")

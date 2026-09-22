# The band is symmetric about 1/2: n_min + n_max = 1, which is just |a1|^2+|a2|^2 = 1. If that
# holds at ARBITRARY contact geometry, then at every fold contact the phase-averaged occupation
# is exactly one half, whatever the metric is doing there. At the bang that is a statement about
# production. At the singular contact A.15 identifies, it would be a statement about TRANSIT:
# how much of what reaches the contact goes through.
#
# Test the symmetry at power-law contacts D(tau) = A sgn(tau)|tau|^n, not just the linear one.

A <- 1
prop <- function(p,n,T,N){ h<-2*T/N
  D<-function(t) A*sign(t)*abs(t)^n
  f<-function(t,M) rbind(-1i*(D(t)*M[1,]+p*M[2,]), -1i*(p*M[1,]-D(t)*M[2,]))
  M<-diag(2)+0i; t<--T
  for(s in seq_len(N)){k1<-f(t,M);k2<-f(t+h/2,M+h/2*k1);k3<-f(t+h/2,M+h/2*k2);k4<-f(t+h,M+h*k3)
    M<-M+h/6*(k1+2*k2+2*k3+k4); t<-t+h}; M }
vecU <- function(t,n,p){ D<-A*sign(t)*abs(t)^n; E<-sqrt(D^2+p^2)
  v<-c(p+0i,-D+E); v/sqrt(sum(Mod(v)^2)) }

cat("   contact exponent n, half-sweep amplitudes from the contact to +T:\n\n")
cat("     n      p      |a1|^2+|a2|^2     n_min      n_max    n_min+n_max   <n> over phase\n")
for (n in c(0.5, 1.0, 2.0)) {
  T <- if (n <= 1) 12 else 6
  for (p in c(0.4, 0.9)) {
    Uh <- prop(p, n, T, 300000)          # full sweep propagator
    Ph <- prop(p, n, T, 300000)
    # half sweep: integrate from 0 to +T by propagating basis columns
    h <- T/150000
    D<-function(t) A*sign(t)*abs(t)^n
    f<-function(t,M) rbind(-1i*(D(t)*M[1,]+p*M[2,]), -1i*(p*M[1,]-D(t)*M[2,]))
    M<-diag(2)+0i; t<-0
    for(s in seq_len(150000)){k1<-f(t,M);k2<-f(t+h/2,M+h/2*k1);k3<-f(t+h/2,M+h/2*k2);k4<-f(t+h,M+h*k3)
      M<-M+h/6*(k1+2*k2+2*k3+k4); t<-t+h}
    uT <- vecU(T,n,p)
    a1 <- sum(Conj(uT)*(M%*%c(1,0))); a2 <- sum(Conj(uT)*(M%*%c(0,1)))
    s2 <- Mod(a1)^2 + Mod(a2)^2
    nmn <- (Mod(a1)-Mod(a2))^2/2; nmx <- (Mod(a1)+Mod(a2))^2/2
    ph  <- seq(0, 2*pi, length.out=4001)[-4001]   # drop the duplicate endpoint: 0 and 2pi are one phase
    avg <- mean(Mod(a1 + exp(1i*ph)*a2)^2/2)
    cat(sprintf("  %4.1f %6.2f %16.9f %10.6f %10.6f %13.9f %15.9f\n",
        n, p, s2, nmn, nmx, nmn+nmx, avg))
  }
}
cat("\n=== reading\n\n")
cat("  n_min + n_max = 1 and the phase average is 1/2 at every exponent tested, to nine\n")
cat("  figures. That is unitarity and nothing else, so it does not care what the geometry at\n")
cat("  the contact is doing. Combined with A.18's result that Theta forces the mode half into\n")
cat("  each component at any contact with an odd diagonal, the statement is:\n\n")
cat("     at every contact the fold makes, the phase-averaged occupation is exactly one half.\n\n")
cat("  At the bang that is about production. At the singular contact of A.15 it would be about\n")
cat("  TRANSIT, and a decohered phase there would give exactly half transmission. STATED AS A\n")
cat("  ROUTE, NOT A RESULT: applying it needs a mode equation for the Kasner interior that this\n")
cat("  paper does not derive, and A.16 separately shows the fold cannot drop the Hawking\n")
cat("  temperature, so nothing here licenses replacing evaporation with transit.\n")

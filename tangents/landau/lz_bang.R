# At a radiation bang a(eta) = a_1 eta, so the conformal Dirac mass term is m a = gamma eta
# and a fermion mode obeys  i d psi/d eta = [[gamma eta, p],[p, -gamma eta]] psi.
# That is the Landau-Zener avoided crossing, exactly, with gap parameter p and sweep rate
# gamma. Nobody in this literature seems to say so. Two payoffs if it holds:
#   (a) LZ gives the occupation in closed form, so the "imported profile" is forced;
#   (b) the FULL sweep needs eta to run from -infinity, which is what the second sheet
#       supplies. A one-sheet universe starts at eta = 0 and gets only half the crossing.
# Test both by direct integration. gamma = 1 throughout.

gamma <- 1
rk4 <- function(p, eta0, eta1, psi0, nstep) {
  h <- (eta1-eta0)/nstep
  f <- function(eta, psi) {
    # -i H psi, H = [[g eta, p],[p, -g eta]]
    rbind(-1i*(gamma*eta*psi[1,] + p*psi[2,]),
          -1i*(p*psi[1,] - gamma*eta*psi[2,]))
  }
  psi <- psi0; eta <- eta0
  for (s in seq_len(nstep)) {
    k1 <- f(eta,       psi)
    k2 <- f(eta+h/2,   psi+h/2*k1)
    k3 <- f(eta+h/2,   psi+h/2*k2)
    k4 <- f(eta+h,     psi+h*k3)
    psi <- psi + h/6*(k1+2*k2+2*k3+k4)
    eta <- eta + h
  }
  psi
}
# adiabatic eigenvectors of H at eta: lower (E<0) and upper (E>0)
lower <- function(eta,p){ E <- sqrt((gamma*eta)^2+p^2); v <- rbind(p+0i, -gamma*eta-E); sweep(v,2,sqrt(colSums(Mod(v)^2)),"/") }
upper <- function(eta,p){ E <- sqrt((gamma*eta)^2+p^2); v <- rbind(p+0i, -gamma*eta+E); sweep(v,2,sqrt(colSums(Mod(v)^2)),"/") }

p  <- c(0.10,0.20,0.30,0.40,0.50,0.60,0.70,0.80,0.90,1.00,1.20,1.40)
Tf <- 18; N <- 260000

cat("=== FULL sweep, eta from -T to +T: the two-sheet crossing\n\n")
psi0 <- lower(-Tf, p)
psiT <- rk4(p, -Tf, Tf, psi0, N)
uT   <- upper(Tf, p)
beta2_full <- Mod(colSums(Conj(uT)*psiT))^2
x2   <- pi*p^2/gamma
cat("      p     |beta|^2 numeric   exp(-pi p^2/g)    ratio      n=(1-sqrt(1-e))/2\n")
for (i in seq_along(p)) {
  cat(sprintf("   %5.2f %18.8f %16.8f %10.5f %20.8f\n",
      p[i], beta2_full[i], exp(-x2[i]), beta2_full[i]/exp(-x2[i]),
      (1-sqrt(1-exp(-x2[i])))/2))
}
cat(sprintf("\n   max relative error against exp(-x^2) over the range: %.3e\n",
    max(abs(beta2_full/exp(-x2)-1))))

cat("\n=== HALF sweep, eta from 0 to +T: the one-sheet start\n\n")
psi0h <- lower(0, p)
psiTh <- rk4(p, 0, Tf, psi0h, N/2)
beta2_half <- Mod(colSums(Conj(upper(Tf,p))*psiTh))^2
cat("      p     |beta|^2 numeric    gamma^2/(16 p^4)    ratio\n")
ph <- c(2.0,2.5,3.0,3.5,4.0)
psi0h2 <- lower(0, ph); psiTh2 <- rk4(ph, 0, Tf, psi0h2, N/2)
b2 <- Mod(colSums(Conj(upper(Tf,ph))*psiTh2))^2
for (i in seq_along(ph))
  cat(sprintf("   %5.2f %18.3e %18.3e %10.4f\n", ph[i], b2[i], gamma^2/(16*ph[i]^4),
      b2[i]/(gamma^2/(16*ph[i]^4))))
ls <- diff(log(b2))/diff(log(ph))
cat(sprintf("\n   log-slope of the half-sweep tail: %s   (p^-4 predicts -4)\n",
    paste(sprintf("%.3f", ls), collapse=", ")))

cat("\n=== which occupation does the abundance use?\n\n")
Ia <- (1/pi^2)*integrate(function(x) x^2*exp(-x^2), 0, Inf)$value
Ib <- (1/pi^2)*integrate(function(x) x^2*(1-sqrt(1-exp(-x^2)))/2, 0, Inf)$value
cat(sprintf("   I from |beta|^2 = exp(-x^2)          : %.7f   (A.1's I_b)\n", Ia))
cat(sprintf("   I from n = (1-sqrt(1-exp(-x^2)))/2   : %.7f   (4.2 and 5.1's I)\n", Ib))
cat(sprintf("   ratio %.4f, and since M_1 ~ I^(-2/5) the mass would move by %.4f,\n",
    Ia/Ib, (Ia/Ib)^(-2/5)))
cat(sprintf("   i.e. 491.6 PeV -> %.1f PeV and the line 245.8 -> %.1f PeV.\n",
    491.6*(Ia/Ib)^(-2/5), 245.8*(Ia/Ib)^(-2/5)))

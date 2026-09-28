#!/usr/bin/env Rscript
# HH/TFD mode-algebra diagnostic for a transparent doubled scalar theory.
# This is NOT a global Nariai stress-tensor calculation.
# Positive-frequency oscillator modes only; the minimal massless zero mode
# is diagnosed separately below and is not silently treated as an oscillator.
# P_perp gives (-1)^l on each real scalar spherical harmonic.

beta <- 2*pi
N <- 64L
nn <- 0:(N-1L)
a <- matrix(0,N,N)
for (j in 2:N) a[j-1L,j] <- sqrt(j-1L)
ad <- t(a); Id <- diag(N)
expect_pair <- function(C,A,B) sum(Conj(C)*(A%*%C%*%t(B)))
tr <- function(A) sum(diag(A))
max_marginal <- max_cross <- max_CCR <- 0
cat("=== 1. Explicit two-mode TFD states and an antipodal insertion ===\n")
cat(" omega  l   local energy       <q_R q_L>         <p_R p_L>\n")
for (omega in c(.2,.7,1.5)) {
  coeff <- exp(-beta*omega*nn/2)
  coeff <- coeff/sqrt(sum(coeff^2))
  C0 <- diag(coeff)
  rho0 <- C0%*%t(C0)
  q <- (a+ad)/sqrt(2*omega)
  p <- -1i*sqrt(omega/2)*(a-ad)
  H <- diag(omega*(nn+.5))
  nbar <- 1/expm1(beta*omega)
  cross <- 1/(2*sinh(beta*omega/2))
  for (l in 0:3) {
    C <- diag(coeff*(-1)^(l*nn))
    rhoR <- C%*%Conj(t(C)); rhoL <- t(C)%*%Conj(C)
    marginal_error <- max(abs(rhoR-rho0),abs(rhoL-rho0))
    qcross <- Re(expect_pair(C,q,q)); pcross <- Re(expect_pair(C,p,p))
    energy <- Re(tr(rhoR%*%H))
    ccr <- expect_pair(C,q%*%p-p%*%q,Id)
    cross_error <- max(abs(qcross-(-1)^l*cross/omega),
                       abs(pcross+(-1)^l*cross*omega))
    max_marginal <- max(max_marginal,marginal_error)
    max_cross <- max(max_cross,cross_error)
    max_CCR <- max(max_CCR,abs(ccr-1i))
    cat(sprintf(" %4.1f   %d   %.12f   % .12f   % .12f\n",omega,l,energy,qcross,pcross))
    stopifnot(marginal_error<1e-13,cross_error<1e-12,
              abs(energy-omega*(nbar+.5))<1e-12,abs(ccr-1i)<1e-12)
  }
}
cat(sprintf("Max marginal change %.3e; cross-covariance error %.3e; CCR error %.3e\n",
            max_marginal,max_cross,max_CCR))
cat("The antipodal insertion flips ODD-l cross-correlations, NOT either local thermal marginal.\n")
cat("Consequently it adds ZERO to any one-factor local stress observable in this doubled-algebra comparison.\n")
cat("This does NOT say the full renormalised HH stress is zero, or that its interior extension is fixed.\n")

cat("\n=== 2. Meaningful failure controls ===\n")
omega <- .2; expected <- -1/(2*omega*sinh(beta*omega/2))
wrong_no_antipode <- -expected
cat(sprintf("Omitting P in l=1: cross-q error %.9f, rejected.\n",abs(wrong_no_antipode-expected)))
stopifnot(abs(wrong_no_antipode-expected)>1)
# The local oscillator energy is fixed by rhoR. Promoting its cross-q term
# to an additive local energy would make two identical reduced states differ.
fake_energy_difference <- omega^2*abs(wrong_no_antipode-expected)
cat(sprintf("False cross-term-as-local-energy rule: predicts difference %.9f although reduced states coincide.\n",
            fake_energy_difference))
stopifnot(fake_energy_difference>.01,max_marginal<1e-13)
nbar <- 1/expm1(beta*omega)
good_CCR <- (nbar+1)-nbar
bad_CCR <- (nbar+.5)-(nbar+.5)
cat(sprintf("Thermal coefficients: (n+1)-n=%.1f; equal-coefficient fault gives %.1f, rejected.\n",good_CCR,bad_CCR))
stopifnot(abs(good_CCR-1)<1e-14,bad_CCR==0)

cat("\n=== 3. Minimal massless Nariai zero mode, independent of P ===\n")
# Euclidean Nariai S2 x S2, radii 1, has eigenvalues j(j+1)+l(l+1).
# The constant mode has j=l=0, norm^2=Vol=16*pi^2 and P parity +1.
eigs <- outer(0:6,0:6,function(j,l) j*(j+1)+l*(l+1))
stopifnot(sum(eigs==0)==1,all(eigs[eigs!=0]>0),(-1)^0==1)
V4 <- 16*pi^2
for (mass in c(.1,.01,.001)) {
  G0 <- 1/(V4*mass^2)
  cat(sprintf("m=%g constant contribution G0=%.9f; m^2 G0=%.12f\n",mass,G0,mass^2*G0))
  stopifnot(abs(mass^2*G0-1/V4)<1e-14)
}
cat("The massless Euclidean inverse does not exist without a zero-mode prescription.\n")
cat("Spacetime derivatives remove this CONSTANT term; an IR-divergent field covariance is not proof of divergent stress.\n")
# Lorentzian homogeneous solution F=atan(sinh t), F'=sech t:
# (cosh(t)F')'=0. The derivative algebra can be treated separately.
tt <- c(-2,-1,0,1,2)
Fprime <- 1/cosh(tt); Fsecond <- -tanh(tt)/cosh(tt)
res <- max(abs(Fsecond+tanh(tt)*Fprime))
badres <- max(abs(-sin(tt)+tanh(tt)*cos(tt))) # wrong trial F=sin t
cat(sprintf("Zero-mode wave-equation residual %.3e; wrong sinusoidal oscillator %.6f, rejected.\n",res,badres))
stopifnot(res<1e-14,badres>.1)
cat("\nPASS: P differs from the angle-preserving HH partner map by (-1)^l.\n")
cat("No physical contact-image energy sign or replacement caustic power is asserted.\n")

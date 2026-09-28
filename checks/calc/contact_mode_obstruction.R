#!/usr/bin/env Rscript
# Ordinary linear Z2 scalar projection under the TIME-REVERSING Nariai fold.
# This tests a quantisation premise, NOT the sign of any physical stress tensor.
# It does not exclude antiunitary gluing or observer-patch quantisation.
# Global metric (radius 1): -dt^2+cosh(t)^2 dtheta^2+dOmega_2^2.
# N:(t,theta,n)->(-t,pi-theta,-n).
# q''+tanh(t)q'+[k^2 sech(t)^2+l(l+1)+m^2+4xi]q=0.
# E(0)=1,E'(0)=0; O(0)=0,O'(0)=1.

evolve <- function(t,k,l,m2=1,xi=0,steps=2400L,damping_sign=1) {
  h <- t/steps; y <- c(1,0,0,1)
  fun <- function(x,v) {
    potential <- k^2/cosh(x)^2+l*(l+1)+m2+4*xi
    c(v[2],-damping_sign*tanh(x)*v[2]-potential*v[1],
      v[4],-damping_sign*tanh(x)*v[4]-potential*v[3])
  }
  for (j in seq_len(steps)) {
    x <- (j-1)*h
    a <- fun(x,y); b <- fun(x+h/2,y+h*a/2)
    c1 <- fun(x+h/2,y+h*b/2); d <- fun(x+h,y+h*c1)
    y <- y+h*(a+2*b+2*c1+d)/6
  }
  y
}

cat("=== 1. Nariai scalar mode transport, parity, and KG norm ===\n")
worst_norm <- 0; worst_parity <- 0
for (k in 0:3) for (l in 0:3) for (t in c(.3,1,pi/2,2)) {
  p <- evolve(t,k,l); n <- evolve(-t,k,l)
  W <- cosh(t)*(p[1]*p[4]-p[2]*p[3])
  norm_error <- abs(W-1)
  parity_error <- max(abs(n-c(p[1],-p[2],-p[3],p[4])))
  worst_norm <- max(worst_norm,norm_error)
  worst_parity <- max(worst_parity,parity_error)
  stopifnot(norm_error<1e-8,parity_error<1e-11)
}
cat(sprintf("64 (k,l,t) cases: max KG-Wronskian drift %.3e; parity error %.3e\n",
            worst_norm,worst_parity))
cat("Cosine spatial parity d=(-1)^(k+l); sine parity d=(-1)^(k+l+1).\n")
cat("N-parities: E*d has d; O*d has -d. Keeping one parity drops its canonical partner.\n")

cat("\n=== 2. Entire finite-dimensional parity sectors, not just single modes ===\n")
# At t=0, N acts on initial data (q,p) as (D q,-D p).
# Any spatial involution D with both parities gives the same algebraic conclusion.
degree <- 8L; D <- diag(c(1,-1,-1,1,-1,1,1,-1))
zero <- matrix(0,degree,degree); ident <- diag(degree)
Omega <- rbind(cbind(zero,ident),cbind(-ident,zero))
S <- rbind(cbind(D,zero),cbind(zero,-D))
stopifnot(max(abs(t(S)%*%Omega%*%S+Omega))==0)
for (eta in c(1,-1)) {
  P <- (diag(2*degree)+eta*S)/2
  restricted <- t(P)%*%Omega%*%P
  cat(sprintf("Time-reversing sector eta=%+d: dimension=%d, restricted symplectic rank=%d\n",
              eta,qr(P)$rank,qr(restricted)$rank))
  stopifnot(qr(P)$rank==degree,max(abs(restricted))==0)
}
Pplus <- (diag(2*degree)+S)/2; Pminus <- (diag(2*degree)-S)/2
stopifnot(qr(t(Pplus)%*%Omega%*%Pminus)$rank==degree)
cat("Opposite parity sectors pair nondegenerately: the full cover has not lost its canonical degrees of freedom.\n")

cat("\n=== 3. A mode-level image-kernel failure ===\n")
# A finite mode test, not a claimed Hadamard state: q=(E-iO)/sqrt(2) has unit KG norm.
# For a k=l=0 mode, N acts as time reflection.
tt <- c(-1,-.4,.4,1)
v <- vapply(tt,function(t) {
  a <- evolve(t,0,0); (a[1]-1i*a[3])/sqrt(2)
},complex(1))
W <- outer(v,Conj(v)); J <- diag(4)[4:1,]
stopifnot(max(abs(W-Conj(t(W))))<1e-12,max(abs(Im(W)))>.1)
for (eta in c(1,-1)) {
  one <- W+eta*W%*%J
  Q <- diag(4)+eta*J
  two <- Q%*%W%*%t(Q)/2
  herm_error <- max(abs(one-Conj(t(one))))
  comm_size <- max(abs(two-t(two)))
  cat(sprintf("eta=%+d one-sided image Hermiticity error %.6f; two-sided commutator size %.3e\n",
              eta,herm_error,comm_size))
  stopifnot(herm_error>.1,comm_size<1e-11,
            min(eigen(two,symmetric=FALSE,only.values=TRUE)$values |> Re())>-1e-10)
}
cat("One-sided image addition can fail Hermiticity; two-sided positive projection loses the CCR.\n")

cat("\n=== 4. Failure controls ===\n")
# A purely SPATIAL antipodal quotient is symplectic and supports canonical parity sectors.
Sspace <- rbind(cbind(D,zero),cbind(zero,D))
stopifnot(max(abs(t(Sspace)%*%Omega%*%Sspace-Omega))==0)
for (eta in c(1,-1)) {
  P <- (diag(2*degree)+eta*Sspace)/2
  rk <- qr(t(P)%*%Omega%*%P)$rank
  cat(sprintf("Spatial-only control eta=%+d: nonzero restricted rank %d\n",eta,rk))
  stopifnot(rk==degree)
}
bad <- evolve(1.5,1,2,damping_sign=-1)
badW <- cosh(1.5)*(bad[1]*bad[4]-bad[2]*bad[3])
cat(sprintf("Wrong damping sign: conserved-Wronskian test returns %.9f, not 1.\n",badW))
stopifnot(abs(badW-1)>1)
cat("\nPASS. Ordinary linear parity-mode quantisation fails for this time-reversing fold.\n")
cat("The physical image stress SIGN is NOT computed or fixed by this result.\n")
cat("Next route: specify local/antiunitary or doubled-algebra state, then apply point-split stress to its global caustic-resolved kernel.\n")

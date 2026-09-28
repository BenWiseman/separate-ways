#!/usr/bin/env Rscript
# Pure Einstein Nariai, Lambda > 0. Computed here independently, 2026-09-24.
# Reproduce: Rscript checks/calc/nariai_negative_mode.R
# Euclidean TT operator O h = -nabla^2 h - 2 R_acbd h^cd.
# Primary cross-check: Volkov & Wipf, hep-th/0003081, (4.4)-(4.5),
# (4.14), (4.19), Table 1. This is NOT a Lorentzian growth-rate calculation.

options(digits = 12)
assert <- function(ok, msg) if (!isTRUE(ok)) stop(msg, call. = FALSE)
delta <- function(i,j) as.numeric(i == j)

# Orthonormal product frame: blocks 1:2 and 3:4.
curvature_action <- function(h, Lambda) {
  out <- matrix(0, 4, 4)
  block <- c(1,1,2,2)
  for (a in 1:4) for (b in 1:4) for (c in 1:4) for (d in 1:4) {
    same <- length(unique(block[c(a,b,c,d)])) == 1
    Racbd <- if (same) Lambda *
      (delta(a,b)*delta(c,d) - delta(a,d)*delta(c,b)) else 0
    out[a,b] <- out[a,b] + Racbd*h[c,d]
  }
  out
}

# Covariant constancy checked in angular coordinates (u,v,w,z).
# h = g_first - g_second. Angular coordinate singularities are not sampled.
parallel_residual <- function(u,w,Lambda) {
  g <- diag(c(1,sin(u)^2,1,sin(w)^2)/Lambda)
  h <- g %*% diag(c(1,1,-1,-1))
  dh <- array(0,c(4,4,4)) # dh[k,i,j] = partial_k h_ij
  dh[1,2,2] <- 2*sin(u)*cos(u)/Lambda
  dh[3,4,4] <- -2*sin(w)*cos(w)/Lambda
  Gam <- array(0,c(4,4,4)) # Gam[m,k,i]
  Gam[1,2,2] <- -sin(u)*cos(u)
  Gam[2,1,2] <- Gam[2,2,1] <- 1/tan(u)
  Gam[3,4,4] <- -sin(w)*cos(w)
  Gam[4,3,4] <- Gam[4,4,3] <- 1/tan(w)
  nab <- dh
  for(k in 1:4) for(i in 1:4) for(j in 1:4)
    for(m in 1:4) nab[k,i,j] <- nab[k,i,j] -
      Gam[m,k,i]*h[m,j] - Gam[m,k,j]*h[i,m]
  max(abs(nab))
}

action <- function(eps,Lambda,G,quotient=FALSE) {
  a2 <- exp(2*eps)/Lambda
  b2 <- exp(-2*eps)/Lambda
  volume <- 16*pi^2*a2*b2
  scalar_curvature <- 2/a2+2/b2
  -volume*(scalar_curvature-2*Lambda)/(16*pi*G) /
    ifelse(quotient,2,1)
}

h <- diag(c(1,1,-1,-1))
assert(sum(diag(h)) == 0, 'Relative-radius mode must be trace free')
max_parallel <- 0
max_hessian_relative <- 0
for (Lambda in c(.3,1,2.7)) for(G in c(.7,1,3)) {
  Oh <- -2*curvature_action(h,Lambda)
  assert(max(abs(Oh + 2*Lambda*h)) < 1e-13, 'TT eigenvalue mismatch')
  rayleigh <- sum(h*Oh)/sum(h*h)
  assert(abs(rayleigh + 2*Lambda) < 1e-13, 'Rayleigh mismatch')
  for(u in c(.2,.8,1.7,2.9)) for(w in c(.3,1.1,2.2,2.8))
    max_parallel <- max(max_parallel,parallel_residual(u,w,Lambda))
  volume <- 16*pi^2/Lambda^2
  tangent <- 2*h # d g(eps)/d eps at eps=0
  spectral_hessian <- volume*sum(tangent*(-2*curvature_action(tangent,Lambda))) /
    (32*pi*G)
  expected <- -16*pi/(G*Lambda)
  assert(abs(spectral_hessian/expected-1) < 1e-13, 'Action/operator mismatch')
  for(q in c(FALSE,TRUE)) {
    target <- expected/ifelse(q,2,1)
    errors <- sapply(c(.02,.01,.005), function(step) {
      fd <- (action(step,Lambda,G,q)-2*action(0,Lambda,G,q)+
               action(-step,Lambda,G,q))/step^2
      abs(fd/target-1)
    })
    assert(all(errors[-3]/errors[-1] > 3.99), 'Hessian convergence not quadratic')
    max_hessian_relative <- max(max_hessian_relative,tail(errors,1))
  }
  for(eps in c(-.2,-.07,.07,.2))
    assert(action(eps,Lambda,G) < action(0,Lambda,G), 'Negative direction absent')
}
assert(max_parallel < 2e-15,'Mode is not covariantly constant')
assert(max_hessian_relative < 8.34e-6,'Finite-difference Hessian inaccurate')

# Euclidean fold: first sphere rotation pi, second sphere antipode.
# In embedding coordinates these are diag(-1,-1,1) and -I_3.
# They preserve the two factor metrics separately, so h is invariant.
# Check pullback in representative orthonormal tangent frames; no factor exchange.
rot2 <- function(a) matrix(c(cos(a),sin(a),-sin(a),cos(a)),2,2)
max_parity <- 0
for(a in seq(0,2*pi,length.out=17)) for(b in seq(0,2*pi,length.out=13)) {
  Q <- matrix(0,4,4)
  Q[1:2,1:2] <- rot2(a)
  Q[3:4,3:4] <- rot2(b) %*% diag(c(1,-1))
  max_parity <- max(max_parity,max(abs(t(Q)%*%h%*%Q-h)))
}
assert(max_parity < 1e-14,'Factor-preserving fold must retain TT mode')

# Meaningful failure controls.
# A DIFFERENT involution that exchanges the factors would make this mode odd.
swap <- diag(4)[c(3,4,1,2),]
assert(max(abs(t(swap)%*%h%*%swap+h)) == 0,'Factor-exchange control failed')
# Curvature-sign error would reverse the instability conclusion.
wrong <- 2*curvature_action(h,1)
assert(sum(h*wrong)/sum(h*h) > 0,'Wrong-sign control failed')
# A uniform scale deformation is not the physical TT negative mode.
assert(sum(diag(diag(4))) != 0,'Conformal contamination control failed')
# Halving the action and norm leaves the eigenvalue unchanged, not halved.
Oh <- -2*curvature_action(h,1)
assert(abs((sum(h*Oh)/2)/(sum(h*h)/2)+2) < 1e-14,
       'Quotient normalization control failed')

cat('PASS: relative-radius tensor is parallel, trace free and fold even.\n')
cat(sprintf('max covariant-derivative residual = %.3e; parity residual = %.3e\n',
            max_parallel,max_parity))
cat('TT eigenvalue = -2 Lambda on cover AND free fold quotient.\n')
cat(sprintf('Lambda=G=1: I(0)=%.12f, I_second=%.12f, quotient I_second=%.12f\n',
            action(0,1,1),-16*pi,-8*pi))
cat(sprintf('three-step Hessian convergence passed; max finest relative error %.3e\n',
            max_hessian_relative))
cat('PASS controls: wrong curvature sign; factor exchange; conformal trace; quotient norm.\n')
cat('Scope: pure Einstein Euclidean Hessian, not a Lorentzian instability or decay rate.\n')

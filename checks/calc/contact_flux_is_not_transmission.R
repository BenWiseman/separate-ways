#!/usr/bin/env Rscript
# Lane 4: contact access is not inter-sheet transmission or mass enhancement.
# Schwarzschild M=1, TEST particles / stationary conserved number flux.
# No metric backreaction or microscopic matter-transfer law is supplied here.
# Conditional on an inward future horizon crossing, not all particles sent
# towards a hole from infinity (which would require a capture distribution).

cat("=== 1. Every captured classical geodesic reaches the contact radius ===\n")
# Normalisation g(u,u)=-eps and Killing constants E,L give
# rdot^2=E^2-f(r)*(eps+L^2/r^2), f=1-2/r.
# For 0<r<2, f<0. An inward timelike ray has no radial turning point.
radii <- seq(1,1.999,length.out=301)
count <- 0; min_radicand <- Inf
for (E in c(.2,.7,1,2)) for (L in c(0,.5,2,10)) for (eps in c(0,1)) {
  rad <- E^2-(1-2/radii)*(eps+L^2/radii^2)
  stopifnot(all(rad>0))
  min_radicand <- min(min_radicand,rad); count <- count+1L
}
cat(sprintf("%d geodesic families x %d radii: all radicands positive, minimum %.9f\n",
            count,length(radii),min_radicand))
cat("Conditional contact-access fraction = 1. The smaller sphere is crossed later, not bypassed by most infall.\n")
# A deliberately false sign for f must make some alleged infall imaginary.
wrong_rad <- .2^2-abs(1-2/radii)*(1+10^2/radii^2)
cat(sprintf("Planted exterior-sign fault: %d negative radicands, rejected.\n",sum(wrong_rad<0)))
stopifnot(any(wrong_rad<0))

cat("\n=== 2. Conserved flux, not an area-fraction probability ===\n")
rr <- c(1.99,1.7,1.3,1)
Jr <- -1/(4*pi*rr^2)
flux <- -4*pi*rr^2*Jr
stopifnot(max(abs(flux-1))<1e-14)
cat("   r      inward current Jr      integrated inward flux\n")
for (j in seq_along(rr)) cat(sprintf(" %.2f   % .12f      %.12f\n",rr[j],Jr[j],flux[j]))
# Keeping flux density constant while shrinking area is NOT conserved infall.
wrong_Jr <- -1/(4*pi*2^2)
wrong_flux_contact <- -4*pi*wrong_Jr
wrong_divergence <- 2*wrong_Jr/rr
cat(sprintf("Area-ratio guess gives %.2f instead of 1; its current has nonzero divergence, max %.9f.\n",
            wrong_flux_contact,max(abs(wrong_divergence))))
stopifnot(abs(wrong_flux_contact-1)>.5,max(abs(wrong_divergence))>.01)

cat("\n=== 3. Geometry supplies an exposure time, not a conversion probability ===\n")
contact_time <- function(E,L) integrate(function(r) {
  denominator <- E^2*r^3+(2-r)*(r^2+L^2)
  ifelse(r==0,0,r^(3/2)/sqrt(denominator))
},0,1,rel.tol=1e-11,subdivisions=1000L)$value
stopifnot(abs(contact_time(1,0)-sqrt(2)/3)<1e-10,
          abs(contact_time(0,0)-(pi/2-1))<1e-10)
for (E in c(.5,1,2)) for (L in c(0,1,3))
  cat(sprintf("E=%3.1f L=%d: proper time r=1 to singularity = %.12f M\n",E,L,contact_time(E,L)))
cat("Even the exposure time depends on the incoming conserved quantities.\n")

cat("\n=== 4. Fold-symmetric, unitary transfer leaves a continuous freedom ===\n")
swap <- matrix(c(0,1,1,0),2,2); Id <- diag(2)
incoming <- c(1,1)
for (theta in c(0,pi/6,pi/4,pi/2)) {
  U <- cos(theta)*Id-1i*sin(theta)*swap
  prob <- Mod(U)^2
  p <- sin(theta)^2
  # Theta=swap*K, K complex conjugation, is an antiunitary involution.
  # Time reversal covariance means Theta U Theta^-1=U^-1.
  stopifnot(max(abs(Conj(t(U))%*%U-Id))<1e-14,
            max(abs(swap%*%Conj(U)%*%swap-Conj(t(U))))<1e-14,
            max(abs(colSums(prob)-1))<1e-14,
            max(abs(prob%*%incoming-incoming))<1e-14)
  cat(sprintf("theta=%.9f, transfer p=%.9f, output fluxes=(%.9f,%.9f)\n",
              theta,p,(prob%*%incoming)[1],(prob%*%incoming)[2]))
}
cat("Same contact geometry, same fold symmetry, canonical unitary evolution: any p from 0 to 1.\n")
cat("This is a family of transport COUNTERMODELS, not a claimed spacetime seam derivation.\n")
cat("For conservative exchange and equal incident supplies, each output stays 1 for every p.\n")
cat("Thus growth factor=1+p is NOT a consequence of transmission alone.\n")
# Planted accounting fault: retain all own input AND add a fraction of the other
# without deducting the reciprocal export. It duplicates flux unless sourced.
p <- .3; wrong_ledger <- Id+p*swap
out <- wrong_ledger%*%incoming
excess <- sum(out)-sum(incoming)
cat(sprintf("Duplicating ledger at p=.3 gives (%.1f,%.1f), total excess %.1f; conservation test rejects it.\n",
            out[1],out[2],excess))
stopifnot(abs(excess-.6)<1e-14,max(abs(colSums(wrong_ledger)-1))>.2)

cat("\nPASS. Geometry fixes access=1 in this Schwarzschild test-particle scope, not transfer p.\n")
cat("A partial transfer prediction needs a transport law and infall distribution.\n")
cat("An exterior mass-growth prediction also needs a consistent energy ledger and backreacted mass definition.\n")

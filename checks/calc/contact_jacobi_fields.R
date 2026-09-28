#!/usr/bin/env Rscript
# Null image-geodesic Jacobi fields at first antipodal contact.
# Derived here rather than quoted; the working is in the header below.
# No packages. No generated files. Schwarzschild M=L=1; Nariai radius=1.
# Jacobi matrices satisfy A(0)=0, dA/dlambda(0)=I.
# Null-screen Van Vleck determinant is Delta=lambda^2/det(A), BEFORE caustic.

rk4 <- function(fun, end, steps) {
  h <- end/steps
  y <- c(0, 0, 1, 0, 1) # affine distance, J1,J1_phi,J2,J2_phi
  out <- matrix(NA_real_, steps+1L, 6L)
  out[1,] <- c(0,y)
  for (i in seq_len(steps)) {
    x <- (i-1)*h
    k1 <- fun(x,y); k2 <- fun(x+h/2,y+h*k1/2)
    k3 <- fun(x+h/2,y+h*k2/2); k4 <- fun(x+h,y+h*k3)
    y <- y+h*(k1+2*k2+2*k3+k4)/6
    out[i+1L,] <- c(i*h,y)
  }
  colnames(out) <- c("phi","lambda","J1","J1_phi","J2","J2_phi")
  out
}
nariai <- function(phi,y) c(1,y[3],-y[2],y[5],y[4])
# Casals et al., Phys. Rev. D 79 (2009) 124043, Eq. (5.43), Sec. V E.
# https://www.barrywardell.net/Research/Papers/Nariai.pdf
# gamma: sphere arc; eta: timelike dS2 distance. On our null ray eta=gamma.
casals_delta <- function(gamma,eta) (gamma/sin(gamma))*(eta/sinh(eta))
schwarzschild <- function(phi,y, sign=1, angular_factor=3) {
  r <- 1+sin(phi); rp <- cos(phi)
  c(r^2, y[3], 2*rp/r*y[3]-sign*angular_factor/r*y[2],
    y[5], 2*rp/r*y[5]+3/r*y[4])
}

cat("=== First-contact null caustic: Nariai and Schwarzschild ===\n")
casals_errors <- wrong_casals_errors <- numeric(0)
for (angle in c(.2,.7,1.5,2.4,pi-.1,pi-.01)) {
  z <- tail(rk4(nariai,angle,1600L),1)
  numerical <- z[,"lambda"]^2/(z[,"J1"]*z[,"J2"])
  target <- casals_delta(angle,angle)
  # Deliberately replace angular focusing by defocusing in the source formula.
  wrong <- (angle/sinh(angle))^2
  casals_errors <- c(casals_errors,abs(numerical/target-1))
  wrong_casals_errors <- c(wrong_casals_errors,abs(numerical/wrong-1))
}
cat(sprintf("Casals Eq. (5.43): six ray endpoints, max relative error %.3e\n",max(casals_errors)))
cat(sprintf("Planted sin(gamma)->sinh(gamma) mismatch: max relative error %.6g, rejected\n",
            max(wrong_casals_errors)))
stopifnot(max(casals_errors)<1e-8,max(wrong_casals_errors)>.5)
for (name in c("Nariai","Schwarzschild")) {
  fun <- if (name=="Nariai") nariai else schwarzschild
  cat("\n",name,"\n",sep="")
  errors <- numeric(3)
  for (j in seq_along(c(400L,800L,1600L))) {
    n <- c(400L,800L,1600L)[j]
    z <- rk4(fun,pi,n)
    exact <- if (name=="Nariai") sin(z[,1]) else (1+sin(z[,1]))*sin(z[,1])
    errors[j] <- max(abs(z[,"J1"]-exact))
    tail <- z[nrow(z),]
    cat(sprintf("n=%4d J1(pi)=% .4e max analytic error=%.4e J2(pi)=%.12f lambda=%.12f\n",
                n,tail["J1"],errors[j],tail["J2"],tail["lambda"]))
    stopifnot(errors[j]<1e-7, all(z[2:n,"J1"]>0), all(z[2:(n+1),"J2"]>0))
    if (name=="Nariai") stopifnot(max(abs(z[,"J2"]-sinh(z[,1])))<1e-7)
    else {
      p <- z[,1]; r <- 1+sin(p)
      exact2 <- 8+cos(p)/r*(-15*p/2+sin(2*p)/4+4*cos(p)-12)
      stopifnot(max(abs(z[,"J2"]-exact2))<2e-7,
                abs(tail["lambda"]-(4+3*pi/2))<1e-10)
    }
  }
  stopifnot(errors[1]/errors[2]>12,errors[2]/errors[3]>10)
  if (name=="Schwarzschild")
    cat(sprintf("Exact other-screen endpoint: 24+15*pi/2=%.12f\n",24+15*pi/2))
  cat("No conjugate point before pi; first angular zero is AT pi.\n")
  cat(" epsilon    det(A)          Delta             epsilon*Delta\n")
  scaled <- numeric(4)
  for (j in seq_along(c(.1,.03,.01,.003))) {
    ep <- c(.1,.03,.01,.003)[j]
    z <- rk4(fun,pi-ep,1600L); last <- z[nrow(z),]
    detA <- last["J1"]*last["J2"]
    Delta <- last["lambda"]^2/detA
    scaled[j] <- ep*Delta
    cat(sprintf(" %.3f   % .9e   % .9e   % .9e\n",ep,detA,Delta,scaled[j]))
    stopifnot(Delta>0)
  }
  stopifnot(abs(scaled[4]/scaled[3]-1)<.04)
  if (name=="Nariai") stopifnot(abs(scaled[4]/(pi^2/sinh(pi))-1)<.01)
  C <- if (name=="Nariai") pi^2/sinh(pi) else (4+3*pi/2)^2/(24+15*pi/2)
  leading <- function(ep) {
    last <- tail(rk4(fun,pi-ep,3200L),1)
    ep*last[,"lambda"]^2/(last[,"J1"]*last[,"J2"])
  }
  richardson <- 2*leading(.0005)-leading(.001)
  cat(sprintf("Leading Delta=C/epsilon: exact C=%.12f; numerical extrapolation=%.12f\n",C,richardson))
  cat(sprintf("Leading sqrt(Delta)=sqrt(C)/sqrt(epsilon): sqrt(C)=%.12f\n",sqrt(C)))
  stopifnot(abs(richardson/C-1)<1e-6,abs(richardson/(C/2)-1)>.9)
  cat("Here epsilon=pi-phi along the ray, NOT radial contact depth s.\n")
}

cat("\n=== Failure controls ===\n")
# If the focusing tidal eigenvalue is given the defocusing sign, no caustic.
bad <- rk4(function(phi,y) schwarzschild(phi,y,sign=-1),pi,1600L)
bad_end <- bad[nrow(bad),"J1"]
cat(sprintf("Wrong tidal sign: J1(pi)=%.9f, fails required zero.\n",bad_end))
stopifnot(abs(bad_end)>1)
# Wrong factor in curvature, while still focusing, must fail the known Killing Jacobi field.
bad2 <- rk4(function(phi,y) schwarzschild(phi,y,angular_factor=2),pi,1600L)
cat(sprintf("Wrong tidal factor 2 instead of 3: J1(pi)=%.9f, fails zero.\n",tail(bad2[,"J1"],1)))
stopifnot(abs(tail(bad2[,"J1"],1))>.1)
# Endpoint is essential: same correct equation at non-antipodal angle is non-caustic.
good <- rk4(schwarzschild,pi/2,800L)
cat(sprintf("Non-antipodal endpoint pi/2: J1=%.9f, J2=%.9f, determinant nonzero.\n",
            tail(good[,"J1"],1),tail(good[,"J2"],1)))
stopifnot(tail(good[,"J1"],1)>1,tail(good[,"J2"],1)>1)
# A finite nonzero Van Vleck limit is rejected by near-endpoint ratio.
a <- tail(rk4(schwarzschild,pi-.01,1600L),1)
b <- tail(rk4(schwarzschild,pi-.001,1600L),1)
da <- a[,"lambda"]^2/(a[,"J1"]*a[,"J2"])
db <- b[,"lambda"]^2/(b[,"J1"]*b[,"J2"])
cat(sprintf("Delta(eps=.001)/Delta(eps=.01)=%.9f, rejects finite-limit hypothesis.\n",db/da))
stopifnot(db/da>9)
cat("\nPASS. Delta is not finite at antipodal contact in either background.\n")
cat("This invalidates the finite-Delta s^-3 argument; it does NOT determine the replacement stress power or sign.\n")

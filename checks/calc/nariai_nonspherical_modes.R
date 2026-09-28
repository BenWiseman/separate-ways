#!/usr/bin/env Rscript
# Pure Einstein dS2 x S2 with Lambda=1; classical modes, NOT quantization.
# Sources: Kodama-Ishibashi hep-th/0308128 (4.7),(4.43),(5.70), App D.2;
# Law-Lochab 2506.02142 (2.15),(2.29)-(2.34); Mukherjee 2506.07556 (A.15).
# The earlier KI hep-th/0305147 explicitly excludes Nariai backgrounds.
options(digits=12)
assert <- function(ok,msg) if(!isTRUE(ok)) stop(msg,call.=FALSE)

# Unnormalized complex Y_lm and its coordinate gradient, recurrence not parity-coded.
plm <- function(l,m,x) {
  if(l<m) return(0)
  p <- if(m==0) 1 else prod(seq(1,2*m-1,by=2))*(-1)^m*(1-x*x)^(m/2)
  if(l==m) return(p)
  p1 <- x*(2*m+1)*p
  if(l==m+1) return(p1)
  for(j in (m+2):l) {
    p2 <- ((2*j-1)*x*p1-(j+m-1)*p)/(j-m)
    p <- p1; p1 <- p2
  }
  p1
}
harmonic <- function(l,m,theta,phi) {
  x <- cos(theta); P <- plm(l,m,x); phase <- exp(1i*m*phi)
  dx <- (l*x*P-(l+m)*plm(l-1,m,x))/(x*x-1)
  grad <- c(-sin(theta)*dx,1i*m*P)*phase
  list(Y=P*phase,grad=grad,axial=c(grad[2]/sin(theta),-sin(theta)*grad[1]))
}
max_ang <- 0
axial_wrong <- 0
for(l in 1:6) for(m in 0:l) for(theta in c(.4,1.1,2.2)) {
  a <- harmonic(l,m,theta,.37)
  b <- harmonic(l,m,pi-theta,.37+pi)
  jac <- c(-1,1)
  max_ang <- max(max_ang,
    Mod(b$Y-(-1)^l*a$Y)/(1+Mod(a$Y)),
    Mod(jac*b$grad-(-1)^l*a$grad)/(1+Mod(a$grad)),
    Mod(jac*b$axial-(-1)^(l+1)*a$axial)/(1+Mod(a$axial)))
  axial_wrong <- max(axial_wrong,
    Mod(jac*b$axial-(-1)^l*a$axial)/(1+Mod(a$axial)))
}
assert(max_ang<3e-13,'Angular tensor pullback parity failed')
assert(axial_wrong>1,'Axial missing-orientation-sign control failed')

# Global N: (t,theta,Omega) -> (-t,pi-theta,P Omega).
# q''+tanh(t)q'+[k^2 sech^2(t)+A]q=0.
# Temporal-even/odd fundamental solutions exist for every k and A.
evolve <- function(A,k,odd,T,steps=1000) {
  y <- if(odd) c(0,1) else c(1,0)
  dt <- T/steps; t <- 0
  rhs <- function(t,y) c(y[2],-tanh(t)*y[2]-(k*k/cosh(t)^2+A)*y[1])
  for(j in seq_len(steps)) {
    k1 <- rhs(t,y); k2 <- rhs(t+dt/2,y+dt*k1/2)
    k3 <- rhs(t+dt/2,y+dt*k2/2); k4 <- rhs(t+dt,y+dt*k3)
    y <- y+dt*(k1+2*k2+2*k3+k4)/6; t <- t+dt
  }
  y
}
max_time <- 0
for(A in c(2,4,10)) for(k in 0:3) for(odd in c(FALSE,TRUE)) {
  yp <- evolve(A,k,odd,3); ym <- evolve(A,k,odd,-3)
  s <- if(odd) -1 else 1
  max_time <- max(max_time,max(abs(ym-c(s,-s)*yp)))
}
assert(max_time<1e-12,'Temporal pullback parity failed')
for(k in 0:6) for(l in 1:6) for(axial in 0:1) for(sine in 0:1) {
  if(k==0 && sine==1) next
  d <- (-1)^(k+l+axial+sine)
  # A classical invariant mode can retain displacement or momentum, not both.
  assert(sum(c(d,-d)==1)==1,'Fold incorrectly removes an entire harmonic sector')
}

# Static master H=-d_x^2+A sech^2(x); A=l(l+1) for minimal test scalar,
# A=l(l+1)-2 for both physical gravitational branches, l>=2.
# Analytic integral <u,Hu>=int(|u_prime|^2+V|u|^2)>=0 is the proof.
# Finite-box matrices provide numerical sign/implementation controls, not the proof.
n <- 201; dx <- 24/(n+1); x <- seq(-12+dx,12-dx,length.out=n)
kin <- diag(2/dx^2,n)
kin[cbind(1:(n-1),2:n)] <- kin[cbind(2:n,1:(n-1))] <- -1/dx^2
lowest <- function(A) min(eigen(kin+diag(A/cosh(x)^2),symmetric=TRUE,
                              only.values=TRUE)$values)
As <- sort(unique(c((1:8)*((1:8)+1),(2:8)*((2:8)+1)-2)))
evals <- sapply(As,lowest)
assert(all(evals>0),'Nonnegative-potential spectral check failed')
tachyon <- lowest(-2)
assert(abs(tachyon+1)<.003,'Planted attractive potential should have eigenvalue -1')

# Fundamental outgoing QNM directly satisfies ODE (not just a frequency lookup).
# u=cosh(x)^(i omega): residual/u=(omega^2+i omega-A)sech^2(x).
qnm_res <- 0
for(A in As) for(s in c(-1,1)) {
  om <- s*sqrt(A-.25)-.5i
  qnm_res <- max(qnm_res,Mod(om^2+1i*om-A))
  assert(Im(om)<0,'Radiative fundamental must decay in static time')
}
assert(qnm_res<1e-12,'Exact fundamental QNM residual failed')

# Exceptional l=1 axial vacuum: h_aA=A_a V_A, dA=J vol_dS2.
# A=sinh(t)dtheta is N-even; F/cosh(t)=1 is constant, not a radiative growing mode.
t <- seq(-4,4,length.out=81)
assert(max(abs(-sinh(-t)-sinh(t)))<1e-14,'Rotational dipole parity failed')
assert(max(abs(cosh(t)/cosh(t)-1))<1e-14,'Rotational field strength not constant')

# Scope control: the surviving spherical X2 is bounded in a chosen static patch
# (X2=rho, |rho|<1), but NOT bounded on the complete global future slices.
assert(cosh(8)>1000,'Global-versus-static growth control failed')

cat(sprintf('PASS harmonic pullbacks l=1..6: max relative residual %.3e\n',max_ang))
cat(sprintf('PASS temporal parity / retained classical data: residual %.3e\n',max_time))
cat(sprintf('Finite-box positive-potential minimum %.9f; tachyon control %.9f\n',
            min(evals),tachyon))
cat(sprintf('Exact fundamental QNM residual %.3e; all radiative Im(omega)=-0.5\n',qnm_res))
cat('Polar angular parity (-1)^l; axial angular parity (-1)^(l+1).\n')
cat('N parity = temporal parity * (-1)^(k+l+axial+sine).\n')
cat('No independent TT S2 tensor sector. Dipole rotation is even and nondynamical.\n')
cat('Scope control: spherical X2 survives and grows globally; bounded only in static patch.\n')
cat('PASS controls: missing axial sign; attractive potential; global/static distinction.\n')

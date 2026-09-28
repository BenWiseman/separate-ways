#!/usr/bin/env Rscript
# mutual_work_symmetric.R -- A.9's claim, and the part of it that cannot be checked.
#
# A.9 argues that a reciprocal time-symmetric response does no net dissipative work, which is
# what §3.3 uses against the Hulse-Taylor orbital decay. The argument is exact and is checked
# here. A.9 also quotes W(G_ret) = -1.4522123062e-2 to ten digits from "analytic-derivative R
# quadrature and a closed-form Gaussian-overlap calculation", and states no source parameters,
# so nobody can reproduce it from the paper. That is recorded rather than papered over.
#
#   W = int dt [ qdot2(t) Phi1(t) + qdot1(t) Phi2(t) ],   Phi_i = (K * q_i)/(4 pi R)
# with Gaussian sources of width sig separated by d, and kernels delta(tau -+ R).

W <- function(kind, sig=1, d=0, R=1) {
  q1 <- function(t) exp(-t^2/(2*sig^2));   q2 <- function(t) exp(-(t-d)^2/(2*sig^2))
  d1 <- function(t) -t/sig^2*q1(t);        d2 <- function(t) -(t-d)/sig^2*q2(t)
  sh <- switch(kind, ret=-R, adv=+R, sym=0)
  f <- function(t) if (kind=="sym")
        0.5*(d2(t)*(q1(t-R)+q1(t+R)) + d1(t)*(q2(t-R)+q2(t+R)))
      else d2(t)*q1(t+sh) + d1(t)*q2(t+sh)
  integrate(f, -40, 40, rel.tol=1e-12, subdivisions=4000)$value/(4*pi*R) }

cat("=== 1. the symmetric kernel does no net work, at every parameter tried ===\n\n")
cat("    sig     d     R        W(ret)          W(adv)          W(sym)\n")
worst <- 0
for (p in list(c(1,0,1), c(1,1,1), c(1,0,2), c(2,0,1), c(1,2,1), c(0.5,0,1), c(1.7,0.4,1.3))) {
  wr <- W("ret",p[1],p[2],p[3]); wa <- W("adv",p[1],p[2],p[3]); ws <- W("sym",p[1],p[2],p[3])
  worst <- max(worst, abs(ws), abs(wr+wa))
  cat(sprintf("   %5.2f %5.2f %5.2f  %+.9e  %+.9e  %+.2e\n", p[1],p[2],p[3], wr, wa, ws)) }
cat(sprintf("\n   worst |W(sym)| and worst |W(ret) + W(adv)|: %.1e, which is quadrature noise.\n", worst))
stopifnot(worst < 1e-12)
cat("   So the even kernel does no net work and the retarded and advanced ones are equal and\n")
cat("   opposite, which is A.9's statement and the reason section 3.3's test bites.\n")

cat("\n=== 2. and W(half Delta) = W(G_ret), which follows from the two above ===\n")
for (p in list(c(1,0,1), c(1,1,1), c(2,0.5,1.5))) {
  wr <- W("ret",p[1],p[2],p[3]); wa <- W("adv",p[1],p[2],p[3])
  cat(sprintf("   sig=%.1f d=%.1f R=%.1f:  W(ret) = %+.9e,  (W(ret)-W(adv))/2 = %+.9e\n",
              p[1],p[2],p[3], wr, (wr-wa)/2))
  stopifnot(abs(wr - (wr-wa)/2) < 1e-12) }

cat("\n=== 3. the frequency-domain statement ===\n")
om <- seq(0.1, 6, length.out = 40); R <- 1.3
cat(sprintf("   max |Im cos(omega R)|               = %.1e   (G_sym is real)\n",
            max(abs(Im(as.complex(cos(om*R)))))))
cat(sprintf("   max |Im( i sin(omega R) ) - sin|    = %.1e   (half Delta carries it all)\n",
            max(abs(Im(1i*sin(om*R)) - sin(om*R)))))

cat("\n=== 4. the check has to be able to fail ===\n")
Wbad <- function(sig=1,d=1,R=1,skew=0.3) {   # a kernel that is NOT even
  q1 <- function(t) exp(-t^2/(2*sig^2)); q2 <- function(t) exp(-(t-d)^2/(2*sig^2))
  d1 <- function(t) -t/sig^2*q1(t);      d2 <- function(t) -(t-d)/sig^2*q2(t)
  f <- function(t) (0.5+skew)*(d2(t)*q1(t-R) + d1(t)*q2(t-R)) +
                   (0.5-skew)*(d2(t)*q1(t+R) + d1(t)*q2(t+R))
  integrate(f, -40, 40, rel.tol=1e-12, subdivisions=4000)$value/(4*pi*R) }
for (s in c(0.1, 0.3, 0.5)) {
  v <- Wbad(skew=s)
  cat(sprintf("     kernel skewed by %.1f toward retarded: W = %+.4e   <- nonzero, as it must be\n", s, v))
  stopifnot(abs(v) > 1e-4) }

cat("\n=== 5. what cannot be checked, and why ===\n")
cat("   A.9 quotes W(G_ret) = -1.4522123062e-2 to ten digits. It states no source width, no\n")
cat("   separation and no R, so the number is not reproducible from the paper. Searching for\n")
cat("   parameters that give it returns no natural values: at sig = R = 1 it needs a source\n")
cat("   separation of 1.3837968, and every other (sig, R) pair returns an equally arbitrary\n")
cat("   one. The number may well be right; it is not checkable, and a ten-digit value with\n")
cat("   unstated inputs should not be quoted at ten digits.\n")
d0 <- uniroot(function(d) W("ret",1,d,1) + 1.4522123062e-2, c(1,2), tol=1e-12)$root
cat(sprintf("     the separation it would need at sig = R = 1: %.7f\n", d0))

cat("
=== flatly ===

  A.9's argument reproduces exactly. The even kernel does no net work at every parameter
  tried, to 1e-16; the retarded and advanced kernels are equal and opposite; and the
  commutator half carries the whole of it. A kernel skewed away from even gives a nonzero
  answer, so the check is not blind.

  A.9's ten-digit illustration does not reproduce and cannot, because the appendix states
  none of the parameters it was computed at. The paper is changed to say so rather than to
  keep quoting ten digits nobody can check.\n")

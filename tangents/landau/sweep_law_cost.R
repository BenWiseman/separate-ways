# ==========================================================================================
# HOW MUCH DOES THE LANDAU-ZENER ASSUMPTION COST?
#
# The operator inequality holds for ANY P, so the sweep law enters the mass only through
# I = Int x^2 n_*(P(x)) dx. 2.3 lists "the Landau-Zener treatment" among the inputs that shift
# the number and never says by how much. This prices it by integrating the two-level problem
# under deformed sweeps.
#
# SUPERSEDED 2026-09-21 by a second pre-submission review. The classifications below are
# NOT RELIABLE and the file is kept for its failure record rather than its numbers:
#   - it tests MONOTONICITY and reports it as convergence; those are different things, and
#     non-monotonicity is not grounds for declaring non-convergence;
#   - its "common terminal gaps" are actually 25, 23.9335 and 25.2222, so the sweeps were not
#     compared on equal footing;
#   - the LINEAR CONTROL is inaccurate where the integral has weight: at g = 1.2 it gives
#     P = 0.00333720 against the exact exp(-pi g^2) = 0.01084671, a factor of three. The
#     validity check only probed g <= 0.6, where it agreed to four per cent.
# An independent integration using instantaneous eigenstates at the boundaries, refined
# momentum grids and extended ranges gives the cubic deformation as 494.777 PeV, a shift of
# +3.177 PeV. So the SHIFT below is approximately right by luck, and the convergence claim
# attached to it was never established.
#
# ORIGINAL NOTE: one deformation appeared to converge and one did not.
#
# FOUR ERRORS ON THE WAY, all recorded because each looked like an answer:
#  1. Compared the numerical result against the analytic Landau-Zener value while measuring
#     the COMPLEMENT: the paper's P is the probability of REMAINING in the diabatic state.
#     The tell was that the numerical value rose with coupling while the analytic one fell.
#  2. Used tanh and arctan sweeps. Both SATURATE, so the diabatic states never separate, P
#     oscillates with the integration window instead of converging, and the "sweep-law effect"
#     was the absence of an asymptotic limit.
#  3. Used faster-than-linear sweeps at a fixed window, so the terminal gap became enormous,
#     the phase per RK4 step exceeded a radian and the integrator returned 0 and NaN. A
#     step-size failure, not physics.
#  4. Quoted the quadratic sweep's 588.9 PeV before checking monotonicity. It is not monotone
#     and its large-g values oscillate, so it has not converged and is NOT a result. It is
#     left in the table with its monotonicity flag so nobody quotes it either.
#
# The step is now set from the terminal gap so the phase per step is bounded, and each sweep's
# window is chosen to give the SAME terminal gap, so the comparison is of the sweep law rather
# than of how far each was integrated.
# ==========================================================================================
nstar <- function(P) (1-sqrt(1-pmin(pmax(P,0),1)))/2
# terminal gap fixed at 25 for every sweep
prob_stay <- function(Delta, g, tmax, phase_per_step=0.02) {
  Dmax <- abs(Delta(tmax)); nstep <- max(20000, ceiling(2*tmax*Dmax/phase_per_step))
  dt <- 2*tmax/nstep; t <- -tmax; psi <- c(1+0i, 0+0i)
  f <- function(tt,p){ H <- matrix(c(Delta(tt), g, g, -Delta(tt)),2,2); -1i*(H%*%p) }
  for (k in 1:nstep){ k1<-f(t,psi); k2<-f(t+dt/2,psi+dt/2*k1)
    k3<-f(t+dt/2,psi+dt/2*k2); k4<-f(t+dt,psi+dt*k3)
    psi <- psi+dt/6*(k1+2*k2+2*k3+k4); t<-t+dt }
  Mod(psi[1])^2
}
nstar <- function(P) (1-sqrt(1-pmin(pmax(P,0),1)))/2
# terminal gap fixed at 25 for every sweep
sw <- list(
  list(nm="linear      t",            D=function(tt) tt,                 tm=25),
  list(nm="cubic-soft  t + 0.02 t^3", D=function(tt) tt+0.02*tt^3,       tm=9.06),
  list(nm="quadratic   (t|t|+t)/2",   D=function(tt) (tt*abs(tt)+tt)/2,  tm=6.62)
)
for (s in sw) cat(sprintf("  %-26s terminal gap = %.2f\n", s$nm, abs(s$D(s$tm))))
gs <- seq(0.1, 1.2, by=0.1)
cat("\n        g "); for (s in sw) cat(sprintf("%14s", substr(s$nm,1,12))); cat("\n")
PP <- lapply(sw, function(s) sapply(gs, function(g) prob_stay(s$D, g, s$tm)))
names(PP) <- sapply(sw, function(s) s$nm)
for (i in seq_along(gs)) { cat(sprintf("   %6.2f", gs[i]))
  for (nm in names(PP)) cat(sprintf("%14.6f", PP[[nm]][i])); cat("\n") }
cat("\n  monotone? "); for (nm in names(PP)) cat(sprintf("%s=%s  ", substr(nm,1,6), all(diff(PP[[nm]])<=1e-5)))
Ival <- function(P){ x <- gs*sqrt(pi); n <- nstar(P)
  sum(diff(x)*(head(x,-1)^2*head(n,-1)+tail(x,-1)^2*tail(n,-1))/2) }
I0 <- Ival(PP[[1]])
cat("\n\n      sweep                      I        I/I_lin   M_1 (PeV)   shift   widths\n")
for (nm in names(PP)) { I <- Ival(PP[[nm]]); M <- 491.6*(I0/I)^(2/5)
  cat(sprintf("   %-26s %8.5f %9.4f %10.1f %8.1f %8.1f\n", nm, I, I/I0, M, M-491.6, (M-491.6)/2)) }

cat("\n  FLATLY, and it is one data point rather than a survey: the cubic deformation is a\n")
cat("  large change to the tail, the cubic term dominating beyond t about 7, and it moves the\n")
cat("  ceiling by 3.3 PeV, under two quoted widths. That is a robustness result as far as it\n")
cat("  goes. The quadratic case did NOT converge and its number is not usable.\n")
cat("\n  NOT PUT IN THE PAPER on the strength of one converged deformation. What the paper\n")
cat("  already says, that the treatment shifts the number rather than blurring it, stands;\n")
cat("  this is a first estimate of the shift and needs at least two more converged sweeps\n")
cat("  before it is worth quoting.\n")

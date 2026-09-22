# ==========================================================================================
# A METHOD FINDING, not a physics result, and it explains an earlier failure.
#
# sweep_law_cost.R deformed the bang's sweep law and reported convergence it had not
# established. A second review found its linear control off by a factor of three at g = 1.2.
# This isolates why, because the reason generalises to anything using the same approach.
#
# THE DIAGNOSIS: there are TWO convergence parameters and only one was tested.
#   - STEP size. Tested by halving. It converges fast and cleanly: at tmax = 28 the norm
#     error falls 2.9e-7 -> 8.9e-12 over four halvings and P is stable to eight digits.
#     AND THE ANSWER IS STILL WRONG, by 4.3 per cent at g = 0.6 and 22.6 at g = 1.2.
#   - WINDOW tmax. This is where the error lives, and it converges as roughly 1/tmax with an
#     oscillating sign:
#         tmax     14      28      56     112
#         error  12.2%    4.3%    3.1%    1.1%
#     Sub-per-cent accuracy needs tmax of order 200, and the earlier work used 6 to 30.
#
# WHY IT FOOLED ME TWICE. Step-halving is the convergence test one reaches for, and it passes
# perfectly while the answer is wrong. And the window error SHRINKS WITH g, so validating at
# small coupling, where the paper's integral has little weight, looks like vindication.
#
# THE FIX, which is what the reviewer used: project onto the INSTANTANEOUS EIGENSTATES at the
# boundaries rather than reading diabatic amplitudes. The slow 1/tmax tail is an artefact of
# asking for an asymptotic label at finite time; the adiabatic projection removes it.
#
# WHAT THIS LEAVES OPEN. 2.2 states the crossing needs "H real, the diagonal odd across the
# contact and the off-diagonal even". Whether REALITY is load-bearing or merely convenient is
# a good question and is NOT answered here: testing a complex coupling with this method would
# repeat the same error. It needs the adiabatic-projection integrator first.
# ==========================================================================================
prop <- function(g, tmax, nstep) {
  dt <- 2*tmax/nstep; t <- -tmax; psi <- c(1+0i, 0+0i)
  f <- function(tt,p){ H <- matrix(c(complex(real=tt), g+0i, g+0i, complex(real=-tt)),2,2)
                       -1i*(H%*%p) }
  for (k in 1:nstep){ k1<-f(t,psi); k2<-f(t+dt/2,psi+dt/2*k1)
    k3<-f(t+dt/2,psi+dt/2*k2); k4<-f(t+dt,psi+dt*k3)
    psi <- psi+dt/6*(k1+2*k2+2*k3+k4); t<-t+dt }
  c(P=Mod(psi[1])^2, nrm=Mod(psi[1])^2+Mod(psi[2])^2)
}
cat("  STEP convergence at fixed window tmax = 28. Exact LZ: exp(-pi g^2).\n\n")
cat("      g    nstep       P          exact      rel err    ||psi||^2-1\n")
for (g in c(0.6,1.2)) for (ns in c(4e4,1.6e5)) {
  r <- prop(g, 28, ns); ex <- exp(-pi*g^2)
  cat(sprintf("   %5.2f %8.0f %11.8f %11.8f %10.2e %12.2e\n", g, ns, r["P"], ex,
              abs(r["P"]-ex)/ex, r["nrm"]-1))
}
cat("\n  Step-converged to 1e-11 in norm, and still 4 to 23 per cent wrong. The step was\n")
cat("  never the problem.\n")
cat("\n  WINDOW convergence at g = 0.6, step held at phase-per-step 5e-3:\n\n")
cat("      tmax       P          exact      rel err\n")
for (tm in c(14,28,56,112)) {
  ns <- ceiling(2*tm*tm/5e-3); P <- prop(0.6, tm, ns)["P"]; ex <- exp(-pi*0.36)
  cat(sprintf("   %7.0f %11.8f %11.8f %10.2e\n", tm, P, ex, abs(P-ex)/ex))
}
cat("\n  Roughly 1/tmax, with the sign oscillating. That is the whole of the earlier failure.\n")

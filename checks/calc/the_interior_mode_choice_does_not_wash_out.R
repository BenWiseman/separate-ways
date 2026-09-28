#!/usr/bin/env Rscript
# Fork 9, a negative with a route attached: the interior mode choice does NOT wash out at large
# l, so the caustic coefficient depends on it, and neither of the two natural choices is the mode.
#
# WHY. the_thermal_factor_is_a_monodromy.R ended by saying the interior mode cannot be chosen at
# the horizon, because a horizon-fixed mode is a monodromy eigenvector and its image correlator is
# then real, which flat space shows is wrong. The obvious hope after that is that the choice does
# not matter where it is used: the caustic divergence comes from the large-l tail of a fully
# coherent sum, so if the two candidate modes converge with l the leading coefficient is the same
# either way and the sign follows without settling the choice. That hope is measured here and it
# is wrong. The difference goes to ONE, not to zero.
#
# THE TWO CANDIDATES. With 2M = 1, mode A is fixed at the HORIZON by u -> e^{-i k rstar}, which is
# interior_modes_at_nonzero_k.R's mode. Mode B is fixed at the SINGULARITY by the non-logarithmic
# Frobenius solution of r^2(r-1)u'' + r u' + [k^2 r^4/(r-1) - l(l+1)r - 1]u = 0, whose indicial
# equation at r = 0 is (s-1)^2 = 0, so u ~ r - l(l+1)r^2 with a log partner. They are compared by
# their LOGARITHMIC derivative at the caustic r = 1/2, which does not care about normalisation.

fail <- 0
note <- function(ok, what) { if (!ok) { cat("   *** FAILED:", what, "\n"); fail <<- fail + 1 } }

fmet  <- function(r) 1 - 1/r
Vpot  <- function(r, l) fmet(r)*(l*(l+1)/r^2 + 1/r^3)
rstar <- function(r) r + log(1 - r)

modeA <- function(k, l, rstop, N = 60000) {          # from the horizon, in rstar
  rs0 <- rstar(1 - 1e-6); h <- (rstar(rstop) - rs0)/N
  st <- c(1 - 1e-6 + 0i, exp(-1i*k*rs0), -1i*k*exp(-1i*k*rs0))
  d <- function(s) { r <- Re(s[1]); c(fmet(r) + 0i, s[3], (Vpot(r, l) - k^2)*s[2]) }
  for (i in seq_len(N)) {
    k1 <- d(st); k2 <- d(st + h*k1/2); k3 <- d(st + h*k2/2); k4 <- d(st + h*k3)
    st <- st + h*(k1 + 2*k2 + 2*k3 + k4)/6
  }
  r <- Re(st[1]); c(u = st[2], dudr = st[3]/fmet(r))
}
modeB <- function(k, l, rstop, r0 = 1e-5, N = 60000) {   # from the singularity, in r
  h <- (rstop - r0)/N; r <- r0
  y <- r0 - l*(l+1)*r0^2 + 0i; yp <- 1 - 2*l*(l+1)*r0 + 0i
  g <- function(r, y, yp)
    c(yp, ((l*(l+1)*r + 1)*y - r*yp - k^2*r^4*y/(r - 1))/(r^2*(r - 1)))
  for (i in seq_len(N)) {
    q1 <- g(r, y, yp); q2 <- g(r + h/2, y + h*q1[1]/2, yp + h*q1[2]/2)
    q3 <- g(r + h/2, y + h*q2[1]/2, yp + h*q2[2]/2); q4 <- g(r + h, y + h*q3[1], yp + h*q3[2])
    y  <- y  + h*(q1[1] + 2*q2[1] + 2*q3[1] + q4[1])/6
    yp <- yp + h*(q1[2] + 2*q2[2] + 2*q3[2] + q4[2])/6; r <- r + h
  }
  c(u = y, dudr = yp)
}

cat("=== 1. the two candidates at the caustic, by logarithmic derivative ===\n")
cat("   Mode A is fixed at the horizon, mode B at the singularity. If the choice washed out,\n")
cat("   the last column would fall with l. It rises to one.\n")
cat("        k     l      u'/u for A            u'/u for B        relative difference\n")
gaps <- c(); k <- 0.6
for (l in c(2, 4, 8, 16, 32)) {
  A <- modeA(k, l, 0.5); B <- modeB(k, l, 0.5)
  la <- A[["dudr"]]/A[["u"]]; lb <- B[["dudr"]]/B[["u"]]
  g <- Mod(la - lb)/Mod(la); gaps <- c(gaps, g)
  cat(sprintf("   %6.2f %5d %9.4f%+9.4fi %9.4f%+9.4fi %18.4e\n",
              k, l, Re(la), Im(la), Re(lb), Im(lb), g))
}
note(all(diff(gaps) > 0), "the difference GROWS with l rather than falling")
note(gaps[length(gaps)] > 0.99, "and reaches one, so the choice dominates the tail the caustic uses")
cat("   So the caustic coefficient depends on the mode choice and cannot be read without it.\n")

cat("\n=== 2. and mode B is not a candidate at all, because it is real ===\n")
cat("   The Klein-Gordon norm on a constant-r slice is W = f(u* u' - u u*'), which vanishes for\n")
cat("   any real solution. Mode B comes out real at every l tried, to machine precision, so it\n")
cat("   carries zero norm and is a standing wave rather than a particle mode.\n")
cat("        l     |Im(u)|/|u| at r = 1/2      |Im(u'/u)|\n")
for (l in c(2, 8, 32)) {
  B <- modeB(k, l, 0.5)
  cat(sprintf("   %6d %22.3e %18.3e\n", l, abs(Im(B[["u"]]))/Mod(B[["u"]]),
              abs(Im(B[["dudr"]]/B[["u"]]))))
  note(abs(Im(B[["dudr"]]/B[["u"]])) < 1e-9, sprintf("mode B is real at l = %d", l))
}
cat("   That is not surprising once seen: the equation has real coefficients and mode B's\n")
cat("   initial data at the singularity is real, so the solution is real all the way out.\n")

cat("\n=== 3. why the flat-space prescription has no counterpart here ===\n")
cat("   In Milne the mode was fixed at LARGE argument, positive frequency in the proper time,\n")
cat("   which is the adiabatic region. The Schwarzschild interior has no such region. The\n")
cat("   Milne radius is R = exp(kappa rstar), and inside, rstar runs from minus infinity at the\n")
cat("   horizon to zero at r = 0, so R covers (0, 1] and stops. The interior is the small-R end\n")
cat("   of Milne and nothing else: it terminates at a curvature singularity before any\n")
cat("   adiabatic region is reached.\n")
for (r in c(1 - 1e-6, 0.9, 0.5, 0.1, 1e-4)) {
  R <- exp(0.5*rstar(r))
  cat(sprintf("   r = %10.6f   rstar = %11.4f   R = exp(rstar/2) = %.6f\n", r, rstar(r), R))
  note(R > 0 && R <= 1 + 1e-12, sprintf("R stays in (0,1] at r = %g", r))
}
cat("   R = 1 exactly at the singularity, since rstar(0) = 0, and R -> 0 at the horizon.\n")
note(abs(exp(0.5*rstar(1e-12)) - 1) < 1e-10, "R is one at the singularity")

cat("\n=== 4. the plant: the comparison has to be able to say 'the same' ===\n")
cat("   Run mode A against itself from a different starting radius. If the logarithmic\n")
cat("   derivative comparison could not report agreement, the result above would be empty.\n")
A1 <- modeA(k, 8, 0.5, N = 60000); A2 <- modeA(k, 8, 0.5, N = 90000)
d <- Mod(A1[["dudr"]]/A1[["u"]] - A2[["dudr"]]/A2[["u"]])/Mod(A1[["dudr"]]/A1[["u"]])
cat(sprintf("   mode A at two step counts: relative difference %.3e\n", d))
note(d < 1e-6, "plant: the comparison reports agreement when there is agreement")

cat("\n=== 5. where this leaves fork 9, with the route and not only the negative ===\n")
cat("   NEGATIVE: the choice does not wash out, so the sign cannot be read off the tail without\n")
cat("   settling it, and neither natural candidate is the mode. Mode B is real and carries no\n")
cat("   norm. Mode A is a monodromy eigenvector and returns a real image correlator, which the\n")
cat("   flat-space contact pair shows is the wrong answer.\n")
cat("   ROUTE: the state is not a single interior solution at all. It is the Hartle-Hawking\n")
cat("   state, fixed by smoothness at the bifurcation surface in Kruskal coordinates, and the\n")
cat("   interior correlator is the continuation of the exterior thermofield double rather than\n")
cat("   a mode with an interior boundary condition. That is the same object the cross weight\n")
cat("   1/sinh(beta w/2) already encodes, and it is what the next brick has to build: the\n")
cat("   exterior pair of solutions, combined as Unruh combines them, continued inside.\n")
cat("   The calibration that closes it is unchanged: Delta^{1/2} -> 3.9004 M s^{-1/2}.\n")

cat("\n=== 6. and the route sharpens to a condition that is already standard ===\n")
cat("   Milne's mode has a second characterisation that does not mention an asymptotic region:\n")
cat("   H^{(2)}_{i k}(z) is the solution RECESSIVE as z goes to -i infinity, since it carries\n")
cat("   e^{-i z} there while H^{(1)} carries e^{+i z}. An analyticity condition transfers where\n")
cat("   an asymptotic one does not, so ask where the same condition lives at a hole. It needs\n")
cat("   R = exp(kappa rstar) large, hence rstar large and POSITIVE, and inside rstar is bounded\n")
cat("   above by zero. Outside it is not: rstar = r + log(r - 1) there and it runs to infinity.\n")
cat("        region      r            rstar        R = exp(rstar/2)\n")
rso <- function(r) r + log(r - 1)
for (rr in c(1 + 1e-6, 1.5, 3, 100)) 
  cat(sprintf("   exterior %10.4f %13.4f %16.4g\n", rr, rso(rr), exp(0.5*rso(rr))))
for (rr in c(1 - 1e-6, 0.5, 1e-4))
  cat(sprintf("   interior %10.4f %13.4f %16.4g\n", rr, rstar(rr), exp(0.5*rstar(rr))))
note(rso(100) > 100, "rstar is unbounded above outside and bounded by zero inside")
note(rstar(1e-9) < 1e-8, "the interior's rstar reaches zero at the singularity and stops")
cat("   So the condition that fixes the mode sits at LARGE r in the EXTERIOR, which is the\n")
cat("   ordinary boundary condition at infinity, and the interior solution is its continuation\n")
cat("   through the horizon. The route is therefore not exotic: solve outside with the standard\n")
cat("   condition, combine the two exterior families the way Unruh does, and continue in. The\n")
cat("   two regions share only the rstar -> -infinity end, which is the horizon, and that is\n")
cat("   where the matching happens.\n")

if (fail == 0) cat("\n   all checks passed\n") else {
  cat(sprintf("\n   %d CHECK(S) FAILED\n", fail)); quit(status = 1)
}

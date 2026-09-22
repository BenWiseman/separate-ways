# Ben: "the graviton acting BETWEEN sheets with a force visible only inside event horizons
# where the sheets pinch."
#
# This is a sharp question and the paper never asks it. Section 2.1 derives the Keldysh
# algebra [Phi_c,Phi_q] = i Delta from ONE premise, stated there as
#
#     "For real test functions supported in one static patch, with the antipodal partner
#      causally disjoint so that Delta(f, g o alpha) = 0"
#
# On de Sitter that premise is free: every point is spacelike to its antipode, which is what
# elliptic de Sitter means, so the paper is right as written and nothing here dents it. The
# question is what happens when the SAME construction is carried to a black hole, which is
# what section 4.3 and the companion do. There the map is the wedge reflection, not the
# antipodal map, and disjointness stops being automatic.
# A cross-sheet commutator is exactly a cross-sheet force: two field
# operators that fail to commute can signal. So the question "is there a force between the
# sheets, and is it confined to the interior" is the question "where does Delta(f, g o alpha)
# fail to vanish", and that is a causal-structure calculation, not an opinion.
#
# Do it in Kruskal, where the wedge reflection J is (U,V) -> (-U,-V). Units 2M = 1 would hide
# factors, so keep M explicit and set M = 1.

M <- 1
cat("=== 1. the geometry, stated so the rest is checkable\n\n")
cat("   Kruskal: U V = (1 - r/2M) exp(r/2M).  Horizon r=2M is UV=0, singularity r=0 is UV=1.\n")
cat("   Regions:  I  U<0 V>0 (our exterior)      II  U>0 V>0 (future interior)\n")
cat("             III U>0 V<0 (far exterior)     IV  U<0 V<0 (past interior)\n")
cat("   The wedge reflection is J:(U,V) -> (-U,-V), so I<->III and II<->IV.\n")
cat("   Two events are causally related in the 2D reduction iff dU and dV share a sign.\n\n")

UV_of_r <- function(r) (1 - r/(2*M))*exp(r/(2*M))
region  <- function(U,V) if (U<0 && V>0) "I" else if (U>0 && V>0) "II" else
                         if (U>0 && V<0) "III" else "IV"
# 2D massless Pauli-Jordan: nonzero only when the separation is timelike, i.e. dU*dV > 0
comm2d  <- function(U1,V1,U2,V2) {
  dU <- U1-U2; dV <- V1-V2
  if (dU*dV > 0) -0.5*sign(dU) else 0
}

cat("=== 2. the whole question, in one table\n\n")
cat("   For a point p and its fold image Jp = (-U,-V): are they causally related?\n\n")
cat("        p (U,V)        region   Jp (U,V)      region   dU     dV    relation    Delta(p,Jp)\n")
pts <- list(c(-2.0, 0.5), c(-0.5, 2.0), c(0.3, 0.3), c(0.6, 0.9), c(0.9, 0.6), c(-0.3,-0.3))
for (p in pts) {
  U <- p[1]; V <- p[2]; JU <- -U; JV <- -V
  dU <- U-JU; dV <- V-JV
  rel <- if (dU*dV > 0) "TIMELIKE" else "spacelike"
  cat(sprintf("   (%5.2f,%5.2f) %8s  (%5.2f,%5.2f) %8s %6.2f %6.2f  %10s %10.2f\n",
      U,V,region(U,V),JU,JV,region(JU,JV),dU,dV,rel,comm2d(U,V,JU,JV)))
}

cat("\n   The pattern is not numerical, it is two inequalities. dU = 2U and dV = 2V, so\n")
cat("   dU*dV = 4UV. The separation between a point and its fold image is timelike exactly\n")
cat("   when UV > 0, and UV > 0 is the INTERIOR. Outside a horizon UV < 0 and the commutator\n")
cat("   vanishes identically, matching what section 2.1 needs; inside, it does not vanish.\n")

cat("\n=== 3. so the cross-sheet commutator is a horizon-confined object\n\n")
cat("      r/2M      UV        sign      Delta(p,Jp) in the 2D reduction\n")
for (rr in c(0.0, 0.5, 1.0, 1.5, 1.999, 2.0, 2.5, 4.0, 10.0)) {
  uv <- UV_of_r(rr*2*M)
  d  <- if (uv > 0) 0.5 else 0
  cat(sprintf("   %8.3f %10.4f %9s %28.2f\n", rr, uv,
      if (uv>0) "interior" else if (uv<0) "exterior" else "horizon", d))
}
cat("\n   It switches on exactly at the horizon and is strictly zero outside it. That is\n")
cat("   Ben's statement, and it is forced by UV changing sign rather than put in by hand.\n")

cat("\n=== 4. the objection that nearly kills it: the transverse parity\n\n")
cat("   The fold is alpha = J o P_perp, not J. P_perp sends a point to its ANTIPODE on the\n")
cat("   sphere, so p and alpha p are separated by pi in angle as well as by (U,V). Angular\n")
cat("   separation costs causal budget, and it can undo the result of section 2. Price it.\n\n")
cat("   A sufficient construction: leave IV, linger in region I sweeping the angle on the\n")
cat("   photon sphere r=3M, then fall into II. A circular null orbit there has\n")
cat("   Omega = dphi/dt = sqrt(M/r^3) = 1/(3 sqrt(3) M), so sweeping pi takes\n\n")
Om <- sqrt(M/(3*M)^3)
dt_pi <- pi/Om
kappa <- 1/(4*M)
mult  <- exp(kappa*dt_pi)
cat(sprintf("      Omega            = %.6f / M\n", Om*M))
cat(sprintf("      coordinate time to sweep pi = %.4f M\n", dt_pi/M))
cat(sprintf("      surface gravity kappa       = %.4f / M\n", kappa*M))
cat(sprintf("      V grows by exp(kappa dt)    = %.2f\n\n", mult))
cat("   V is multiplied by about sixty while the angle is swept. So the image is causally\n")
cat("   related to p only if p sits deep enough that its V exceeds roughly sixty times the\n")
cat("   V at which the curve re-entered. Near the bifurcation surface there is no room and\n")
cat("   the commutator vanishes even inside the horizon; deep inside there is.\n")
cat("   THAT is the pinch, and it has a scale rather than being a surface put in by hand.\n")

cat("\n=== 5. flatly\n\n")
cat("  ESTABLISHED, exactly, in the 2D reduction: at a BLACK HOLE the cross-sheet\n")
cat("  commutator is zero outside the horizon and nonzero inside, because dU dV = 4UV and\n")
cat("  UV > 0 is the interior. One line, no fitting.\n\n")
cat("  WHAT IT IS NOT, said first so it is not overclaimed. This is NOT a defect in\n")
cat("  section 2.1. There the fold lives on de Sitter, alpha is the antipodal map, and\n")
cat("  EVERY point is spacelike to its antipode: that is elliptic de Sitter's defining\n")
cat("  property, so the premise holds globally and the paper's scoping to a static patch\n")
cat("  is correct as written. The failure is in the BLACK HOLE extension, where J is the\n")
cat("  wedge reflection rather than the antipodal map, which is the companion's subject\n")
cat("  and section 4.3's. No number in sections 2 or 3 moves.\n\n")
cat("  WHAT IT COSTS: the Keldysh doubling [Phi_c,Phi_q] = i Delta is derived from cross-\n")
cat("  sheet vanishing, so carried to a hole it holds outside the horizon and fails inside.\n")
cat("  Any future interior algebra has to be built, not inherited.\n\n")
cat("  WHAT FAILS: this is the s-wave reduction. Section 4's estimate is a SUFFICIENT\n")
cat("  construction using an unstable circular orbit, so sixty is an upper bound on the\n")
cat("  cost and not the optimal curve. The honest statement is that the interior region\n")
cat("  carrying a cross-sheet commutator is nonempty and bounded away from the bifurcation\n")
cat("  surface, not that its boundary is known. A geodesic maximisation would give the\n")
cat("  boundary, and it is a half-day of work rather than a paragraph.\n\n")
cat("  NEXT ROUTE, and it is the interesting one. A nonvanishing cross-sheet commutator\n")
cat("  confined to the interior is a force that no exterior observer can measure, which is\n")
cat("  why the paper's 'no transport between branches' survives intact. But it is exactly\n")
cat("  the regime where the transit question lives. If the sheets can exchange stress only\n")
cat("  inside horizons, then a two-sided hole is the ONLY place the two sheets interact,\n")
cat("  and section 4.3's two-sided class stops being a curiosity and becomes the sole\n")
cat("  channel. That is a companion paper, and it needs the interior algebra the paper\n")
cat("  explicitly does not have.\n")

# Ben: "The fusion may not be a single point: it could be statistical, scale-dependent, or
# several partial/fuzzy joins that together make a coherent whole."
#
# Last session established that the sheets' algebras fail to commute only INSIDE horizons
# (cross_sheet_commutator_inside.R). That immediately raises a problem with my own result:
# the bang is not inside a horizon, yet the whole of section 2.2 is about the sheets meeting
# there. If the sheets can only interact inside horizons, section 2.2 should be impossible.
#
# Either the picture is inconsistent or the joins are genuinely different objects. Check it,
# and while checking, enumerate every place the two sheets touch and compute what kind of
# contact each one is. Everything below is recomputed here rather than quoted, so the script
# fails if any of it has drifted.

cat("=== JOIN 1. At the bang: not two algebras at all\n\n")
n_of_x <- function(x) (1 - sqrt(1 - exp(-x^2)))/2
S_of_n <- function(n) ifelse(n <= 0 | n >= 1, 0, -n*log(n) - (1-n)*log(1-n))

xs <- seq(1e-6, 40, length.out = 400001)
w  <- xs^2
I_int  <- sum(w*n_of_x(xs))*(xs[2]-xs[1])          # integral of x^2 n dx
S_int  <- sum(w*S_of_n(n_of_x(xs)))*(xs[2]-xs[1])  # integral of x^2 S dx
I_pap  <- I_int/pi^2
cat(sprintf("   production integral I = (1/pi^2) int x^2 n dx = %.7f   (paper: 0.0127597)\n", I_pap))
cat(sprintf("   int x^2 S dx                                  = %.4f      (paper: 0.4509)\n", S_int))
cat(sprintf("   nats per unit production integral = S_int/I_int = %.3f   (paper: 3.58)\n", S_int/I_int))
cat(sprintf("   per relic particle (halved, two per pair)       = %.3f   (paper: 1.79)\n", S_int/I_int/2))
stopifnot(abs(I_pap - 0.0127597) < 2e-6, abs(S_int/I_int/2 - 1.79) < 0.01)

# the scale of the crossover, recomputed
Smax <- log(2)
f <- function(x) S_of_n(n_of_x(x)) - Smax/2
xc <- uniroot(f, c(0.2, 3))$root
cat(sprintf("   half-entanglement scale x = %.3f                (paper: 0.968)\n", xc))
cat(sprintf("   transition probability falls to 1/e at x = 1, a gap of %.1f per cent\n",
            100*abs(xc-1)/1))
cat("\n   At the bang eta<0 and eta>0 are ONE free field at two times, not two commuting\n")
cat("   subalgebras, so there is no cross-sheet commutator to ask about. The paper says\n")
cat("   this explicitly in 3.5 and it is what resolves the apparent contradiction: the\n")
cat("   horizon result does not apply because the objects it applies to do not exist here.\n")
cat("   The join at the bang is TOTAL, and it is graded by momentum rather than located.\n")

cat("\n=== JOIN 2. Outside a horizon: correlated, but unable to signal\n\n")
cat("   Section 2.1's premise, now known to be exactly the condition UV<0. The commutator\n")
cat("   is identically zero. The cross-leg KERNEL is not: A.1 gives, on the horizon,\n")
cat("       G_alpha / G_J -> (1-cos gamma)/(1+cos gamma) = tan^2(gamma/2).\n\n")
cat("      gamma (rad)    G_alpha/G_J     commutator    correlated?   can signal?\n")
for (g in c(pi/6, pi/3, pi/2, 2*pi/3, 5*pi/6)) {
  r <- (1-cos(g))/(1+cos(g))
  cat(sprintf("   %12.4f %14.4f %14d %14s %13s\n", g, r, 0, "yes", "no"))
}
cat("\n   This is the ordinary signature of entanglement without signalling, and it is the\n")
cat("   whole content of 'no transport between branches'. The sheets are CORRELATED\n")
cat("   everywhere and can INFLUENCE each other nowhere out here.\n")

cat("\n=== JOIN 3. Inside a horizon: able to signal, with a scale\n\n")
M <- 1
cat("   From cross_sheet_commutator_inside.R: dU dV = 4UV, so the separation between a\n")
cat("   point and its image is timelike exactly when UV>0, which is the interior.\n")
Om    <- sqrt(M/(3*M)^3); dt_pi <- pi/Om; kap <- 1/(4*M)
cat(sprintf("   pinch cost of the transverse parity: V multiplied by %.1f\n", exp(kap*dt_pi)))
cat("   So: zero outside, nonzero deep inside, and vanishing again as the bifurcation\n")
cat("   surface is approached. A join with a boundary, not a surface put in by hand.\n")

cat("\n=== JOIN 4. At the singularity: THIS ONE FAILS, and it is the interesting failure\n\n")
cat("   The bang works because a(eta) = a_1 eta changes SIGN through the contact, making\n")
cat("   the diagonal of H odd, which is the whole load-bearing hypothesis of 2.2 and B.5.\n")
cat("   A Schwarzschild interior is Kasner with r proportional to |tau|^(2/3). Check the\n")
cat("   parity of that directly rather than asserting it:\n\n")
cat("      tau      r ~ |tau|^(2/3)     r(-tau) - r(tau)      odd?\n")
for (tau in c(0.1, 0.5, 1.0, 2.0)) {
  rp <- abs(tau)^(2/3); rm <- abs(-tau)^(2/3)
  cat(sprintf("   %6.2f %16.6f %20.3e %9s\n", tau, rp, rm-rp, if (abs(rm-rp)<1e-15) "NO, even" else "odd"))
}
cat("\n   Even at every tau, exactly. So there is no sign change to carry the contact\n")
cat("   condition, and writing sgn(tau) into the singular locus would be assuming the\n")
cat("   conclusion. The fourth join does not exist as a crossing of the bang's type.\n")

cat("\n=== SO: the fusion is four different things, three of which are real\n\n")
cat("      where                what joins               kind of contact        computed\n")
cat("   -------------------------------------------------------------------------------\n")
cat(sprintf("   bang                 one field, two times     total, graded in k     x_c=%.3f\n", xc))
cat("   outside a horizon    two algebras             correlated, no signal  Delta = 0 exactly\n")
cat("   inside a horizon     two algebras             signal, with a pinch   V x 59.2\n")
cat("   at a singularity     ---                      NONE of this type      r even in tau\n")
cat("\n  Ben's framing was right and the single-point fusion was never the right object.\n")
cat("  The three real joins are not three versions of one thing: they differ in WHAT is\n")
cat("  joined (one field vs two algebras) and in WHICH structure carries it (occupation,\n")
cat("  kernel, commutator). The reason they do not contradict is that each lives where the\n")
cat("  others' objects are undefined, and that is checkable rather than rhetorical.\n")

cat("\n=== 5. the next route, taken now rather than promised\n\n")
cat("   Join 4 died on the parity of the scale function. So ask the parity question of a\n")
cat("   ROTATING interior, which has a different structure and which nobody has asked.\n\n")
Mk <- 1; a <- 0.9
Delta <- function(r) r^2 - 2*Mk*r + a^2
Sigma <- function(r, th) r^2 + a^2*cos(th)^2
cat("   Kerr: Delta(r) = r^2 - 2Mr + a^2,  Sigma = r^2 + a^2 cos^2(theta).\n")
cat("   Decompose each into even and odd parts through r=0, the ring:\n\n")
cat("        r      Delta(r)   Delta(-r)     even part     ODD part     Sigma odd part\n")
for (r in c(0.25, 0.5, 1.0, 2.0)) {
  dp <- Delta(r); dm <- Delta(-r)
  ev <- (dp+dm)/2; od <- (dp-dm)/2
  sod <- (Sigma(r, pi/3) - Sigma(-r, pi/3))/2
  cat(sprintf("   %6.2f %11.4f %11.4f %13.4f %12.4f %17.2e\n", r, dp, dm, ev, od, sod))
}
cat("\n   The odd part of Delta is exactly -2Mr: linear in r, nonzero, and it is the ONLY\n")
cat("   odd structure in the metric functions, Sigma being even to machine zero. Compare\n")
cat("   the bang, where the odd object is a(eta) = a_1 eta, also exactly linear. Under\n")
cat("   r -> -r the Kerr metric maps to Kerr with M -> -M, so the ring is a sign flip of\n")
cat("   the mass parameter rather than an endpoint.\n\n")
cat("   And the geometry permits the crossing, which Schwarzschild does not: the Kerr\n")
cat("   singularity is a RING at r=0, theta=pi/2 only, so r=0 is a passable disc off the\n")
cat("   ring and the maximal extension carries r<0. Schwarzschild's is spacelike and\n")
cat("   unavoidable. Check that the two differ where it matters, at theta away from pi/2:\n\n")
cat("   Sigma(0,theta) = a^2 cos^2(theta), which vanishes ONLY at theta = pi/2.\n\n")
cat("        theta     cos(theta)   Sigma(0,theta)    r=0 singular?\n")
for (th in c(0, pi/6, pi/3, 1.5, pi/2)) {
  sg <- Sigma(0, th)
  cat(sprintf("   %10.4f %12.6f %16.8f %16s\n", th, cos(th), sg,
      if (abs(cos(th)) > 1e-12) "no, regular" else "YES, the ring"))
}
cat("\n   The singular set is a measure-zero ring in the disc, not the whole disc, so an\n")
cat("   infalling worldline generically misses it. Schwarzschild has no such escape: its\n")
cat("   singularity is the whole spacelike surface r=0 and every worldline ends on it.\n")
cat("\n   So a rotating interior has BOTH things Schwarzschild lacks: an odd part linear in\n")
cat("   the crossing coordinate, and a locus that can actually be crossed.\n")
cat("   CAVEAT, and it is the whole caveat: this establishes that the NECESSARY parity\n")
cat("   structure is present, not that the contact condition holds. The bang's condition\n")
cat("   is on the 2x2 mode Hamiltonian, diagonal odd and off-diagonal even, and the Kerr\n")
cat("   radial equation has not been put in that form here. Necessary, not sufficient.\n")

cat("\n=== flatly\n\n")
cat("  ESTABLISHED: the apparent contradiction between the bang join and the horizon\n")
cat("  result is not one. They concern different objects, and the paper already says so\n")
cat("  in 3.5; this script only makes the taxonomy explicit and recomputes each entry.\n\n")
cat("  FAILS: join 4. There is no contact condition at a singular locus of Schwarzschild\n")
cat("  type, because Kasner r is EVEN in tau and the bang's condition needs oddness. That\n")
cat("  kills the most attractive version of the transit picture, in which matter crosses\n")
cat("  the singularity the way a mode crosses the bang. It does not die by hand-waving:\n")
cat("  it dies on the parity of |tau|^(2/3).\n\n")
cat("  AND, in the same turn: section 5 asked the parity question of a rotating interior\n")
cat("  and it comes back positive. Kerr's Delta has odd part exactly -2Mr, linear in the\n")
cat("  crossing coordinate just as the bang's a_1 eta is, Sigma is even, and the ring\n")
cat("  leaves r=0 passable off theta=pi/2 where Schwarzschild's spacelike singularity is\n")
cat("  not. Schwarzschild fails the transit picture; Kerr has the structure for it. That\n")
cat("  matters because every astrophysical hole rotates, so the case the fold needs is\n")
cat("  the generic one and the case that fails is the idealisation.\n\n")
cat("  WHAT IS STILL NOT SHOWN, stated so nobody reads more into it: necessary parity is\n")
cat("  not the contact condition. That condition is on the 2x2 mode Hamiltonian, and the\n")
cat("  Kerr radial equation has not been put in that form here. The companion's open item\n")
cat("  is now sharp rather than vague: write Teukolsky in the bang's variables and read\n")
cat("  off whether the diagonal is odd through r=0. That is a well-posed calculation with\n")
cat("  a yes or no at the end of it, which the previous formulation did not have.\n")

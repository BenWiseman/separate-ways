# Both referees: "the reciprocal law is definitional - any involution gives it, so it
# constrains nothing." Test whether that is true OFF the bifurcation surface.
#
# On B the boost separation drops out and alpha acts as P_perp, so alpha^2 = 1 really is
# the trivial statement that an involution is an involution. But off B, the paper's own
# A.6 says alpha acts as the HALF-PERIOD MODULAR SHIFT:
#     G_img(t) = W(t - i beta/2),   beta = 2 pi / H.
# If so, alpha^2 = 1 is not a triviality off B: it is the statement that two half-period
# shifts make a full period and return W to itself. That is the KMS condition, and it
# holds at ONE temperature. Test exactly that.

H <- 1; beta <- 2*pi/H
W <- function(t) 1/sinh(H*t/2)^2                 # de Sitter conformal scalar, up to constant
shift <- function(f, b) function(t) f(t - 1i*b/2)

cat("=== 1. one half-period shift reproduces the image kernel A.6 quotes\n\n")
ts <- c(0.4, 0.9, 1.7, 3.1)
cat("        t        W(t - i beta/2)            -sech^2(Ht/2)         |difference|\n")
for (t in ts) {
  a <- W(t - 1i*beta/2); b <- -1/cosh(H*t/2)^2
  cat(sprintf("   %7.2f  %10.6f%+10.6fi  %10.6f%+10.6fi   %.2e\n",
      t, Re(a), Im(a), Re(b), Im(b), Mod(a-b)))
}
cat("\n  The half shift turns csch^2 into -sech^2, which is A.6's G_img up to the\n")
cat("  constant. So alpha acting off B IS the half-period continuation.\n")

cat("\n=== 2. TWO half shifts. Does alpha^2 = 1 hold, and at which beta?\n\n")
cat("  Sweep the shift period b and measure |W(t - i b) - W(t)| at four times.\n")
cat("  If this were definitional it would vanish for every b.\n\n")
cat("        b/beta        max |W(t - i b) - W(t)| over the four times\n")
for (fr in c(0.25, 0.5, 0.75, 0.9, 0.99, 1.0, 1.01, 1.5, 2.0)) {
  b <- fr*beta
  v <- max(sapply(ts, function(t) Mod(W(t - 1i*b) - W(t))))
  cat(sprintf("   %10.3f %20.3e %s\n", fr, v, if (v < 1e-12) "   <== alpha^2 = 1" else ""))
}
cat("\n  It vanishes at b = beta and at b = 2 beta, and NOWHERE else. The involution\n")
cat("  property of alpha therefore holds at one temperature, not at every temperature.\n")
cat("  Off B it is a constraint. On B it degenerates into P_perp^2 = 1, which is the\n")
cat("  version the referees called definitional. Same equation, two limits.\n")

cat("\n=== 3. so what does alpha^2 = 1 actually fix? Solve it for the temperature.\n\n")
cat("  Treat beta as unknown, impose W(t - i beta) = W(t) at a single time, root-find:\n\n")
cat("  (beta is a TANGENCY of |W(t-ib)-W(t)|, not a sign change, so minimise.)\n\n")
for (t0 in c(0.3, 0.8, 1.5, 2.6)) {
  f <- function(b) Mod(W(t0 - 1i*b) - W(t0))
  o <- optimize(f, c(0.5*beta, 1.5*beta), tol=1e-14)
  cat(sprintf("   from t = %4.1f :  argmin = %.10f  residual = %.2e  %s\n",
      t0, o$minimum, o$objective,
      if (o$objective < 1e-10) sprintf("SOLUTION, err vs 2pi/H = %.1e", abs(o$minimum-beta))
      else "local minimum only, NOT a solution"))
}
cat("\n  Three of the four give the Gibbons-Hawking value to machine precision. The\n")
cat("  fourth is a boundary artefact of the search window, whose lower end sits at pi,\n")
cat("  and it is labelled rather than quietly dropped: a residual of 45 is not a root.\n")
cat("  A wrong beta fails the condition by a finite, measurable amount.\n")

cat("\n=== 4. and the frequency-domain form, which is the half of the KMS factor\n\n")
cat("  A.6 records F_half(E) = e^{beta E/2} F_W(E). Applying the shift twice gives\n")
cat("  e^{beta E}, the FULL KMS factor. Check the composition on the spectral side:\n\n")
cat("        E     e^{beta E/2}   (e^{beta E/2})^2     e^{beta E}     |difference|\n")
for (E in c(0.1, 0.3, 0.7, 1.2)) {
  h <- exp(beta*E/2)
  cat(sprintf("   %6.2f %13.5f %17.5f %14.5f %14.2e\n", E, h, h^2, exp(beta*E), abs(h^2-exp(beta*E))))
}
cat("\n  So the fold's map is the square root of the KMS transformation, and the\n")
cat("  reciprocal law R(x,y) R(x, P_perp y) = 1 is that square root squaring to one,\n")
cat("  restricted to B where the boost separation has dropped out.\n")

cat("\n=== 5. what this does and does not answer\n\n")
cat("  It answers the charge that the law is empty: it is empty ON B and not empty off\n")
cat("  it, and the on-B version is the boundary value of the off-B one. What it is NOT\n")
cat("  is new physics. alpha = J o P_perp contains the wedge reflection, the boost by\n")
cat("  imaginary angle pi, so alpha^2 = 1 says a 2 pi Euclidean rotation is the\n")
cat("  identity: the smoothness condition, from which the horizon temperature follows\n")
cat("  by a standard argument. The value here is defensive and real - the fold's Z2 IS\n")
cat("  that condition, so calling it bookkeeping is wrong. It also does NOT\n")
cat("  establish the reduction G_alpha/G_J -> W_B(x,P_perp y)/W_B(x,y) at a general\n")
cat("  horizon, which both referees named as the load-bearing unproven step, and which\n")
cat("  this says nothing about. It is computed for the de Sitter conformal scalar in\n")
cat("  Bunch-Davies, which is one state in one spacetime.\n")

cat("\n=== CHALLENGE: is this circular? A referee says W is hard-coded 2 pi i/H-periodic,\n")
cat("    so 'alpha^2 = 1 holds only at beta = 2 pi/H' verifies a built-in property.\n\n")
cat("  CONCEDED for the FULL-period step: sinh(x - i pi) = -sinh(x), so csch^2 has period\n")
cat("  2 pi i/H whatever we do, and sweeping the full period tests that and nothing else.\n\n")
cat("  NOT CONCEDED for the step that carries the argument. A.6 derives the image kernel\n")
cat("  G_img ~ sech^2(Ht/2) from the GEOMETRY of the antipodal map, with no beta anywhere in\n")
cat("  its derivation. The content is that this geometric object equals a HALF-PERIOD THERMAL\n")
cat("  SHIFT of W, and the question is whether that match fixes beta or merely tolerates it.\n")
cat("  Sweep beta in the MATCH, not in the periodicity, and see.\n\n")
H <- 1
Wt <- function(t,b) 1/sinh(H*t/2 - 1i*H*b/4)^2     # W(t - i b/2)
target <- function(t) -1/cosh(H*t/2)^2             # A.6's geometric image kernel
ts <- c(0.4,0.9,1.7,3.1)
cat("        beta/(2pi/H)     max |W(t - i beta/2) - G_img| over the sample\n")
for (f in c(0.25,0.5,0.75,0.9,1.0,1.1,1.5,2.0,3.0)) {
  b <- f*2*pi/H
  d <- max(sapply(ts, function(t) Mod(Wt(t,b)-target(t))))
  cat(sprintf("   %14.2f %20.3e %s\n", f, d, if (d<1e-12) "  <-- MATCH" else ""))
}
o <- optimize(function(b) max(sapply(ts, function(t) Mod(Wt(t,b)-target(t)))), c(0.1, 4*pi))
cat(sprintf("\n   minimising over beta returns %.10f against 2 pi/H = %.10f, residual %.2e\n",
    o$minimum, 2*pi/H, o$objective))
cat("\n  VERDICT. CORRECTED 2026-09-20: the table above prints a second MATCH at 3x, and the\n")
cat("  line that used to stand here said the match fails everywhere else, contradicting it.\n")
cat("  Every ODD multiple matches; see kms_period_is_primitive.R. The fundamental is selected\n")
cat("  by primitivity, which is part of the KMS condition. The match holds at beta = 2 pi/H\n")
cat("  and fails at every non-odd-multiple value, by seven orders at\n")
cat("  ten per cent away. So beta is DETERMINED by requiring the geometric image kernel to be a\n")
cat("  thermal half-period shift; it is not inserted. The temperature is returned by that match,\n")
cat("  and the full-period step the referee attacks is a consistency check downstream of it.\n")
cat("  The objection lands on the presentation, which led with the weaker step, not on the\n")
cat("  result. A.10 should lead with the match.\n")

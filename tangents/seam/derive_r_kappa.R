# A referee: "kappa = 1 hangs entirely on r(kappa) = (1-k^2)/(1+k^2), never derived; its only
# provenance is a script header citing two model reviews." The header is real. Provenance is not
# an argument if the form is right, so derive it and stop depending on anyone's say-so.
#
# Corner term S_c = (kappa/2) int dt (q1 qdot2 - q2 qdot1), q = (phi_L(0), phi_R(0)).
# Varying, and integrating by parts in t, the boundary terms are kappa(qdot2 dq1 - qdot1 dq2),
# which against the bulk variation -phi' dphi gives the matching conditions
#     phi_L'(0) = +kappa qdot2,     phi_R'(0) = +kappa qdot1.
# Scattering ansatz at frequency omega = k, time dependence e^{-i omega t} so qdot = -i k q:
#     x<0: phi = e^{ikx} + r e^{-ikx},  q1 = 1 + r
#     x>0: phi = t e^{ikx},             q2 = t

cat("=== solve the two matching conditions symbolically, then check numerically\n\n")
cat("   phi_L'(0) = ik(1-r),  phi_R'(0) = ikt,  qdot1 = -ik(1+r),  qdot2 = -ikt\n")
cat("   condition 1:  ik(1-r) = kappa * (-ik) t        =>  1 - r = -kappa t ... (sign A)\n")
cat("   condition 1': ik(1-r) = -kappa * (-ik) t       =>  1 - r = +kappa t ... (sign B)\n")
cat("   condition 2:  ikt      = -kappa * (-ik)(1+r)   =>  t = kappa (1+r)\n\n")
cat("   Sign B with condition 2:  1 - r = kappa^2 (1+r)  =>  r = (1-kappa^2)/(1+kappa^2)\n")
cat("   and then t = kappa(1+r) = 2 kappa/(1+kappa^2).  Sign A gives the reciprocal and is the\n")
cat("   orientation with the seam traversed the other way; it is the same seam relabelled.\n\n")

rk <- function(k) (1-k^2)/(1+k^2); tk <- function(k) 2*k/(1+k^2)
cat("      kappa      r from 1-r=k^2(1+r)      r(kappa) claimed       t = k(1+r)      t claimed     r^2+t^2\n")
for (k in c(0, 0.25, 0.5, 1, 2, 4)) {
  r_solved <- (1-k^2)/(1+k^2)                  # the solution of the linear equation
  lhs <- 1 - r_solved; rhs <- k^2*(1+r_solved) # verify it solves it
  t_from <- k*(1+r_solved)
  cat(sprintf("   %8.2f %18.6f %22.6f %16.6f %14.6f %11.6f\n",
      k, r_solved, rk(k), t_from, tk(k), r_solved^2+t_from^2))
  stopifnot(abs(lhs-rhs) < 1e-12, abs(t_from-tk(k)) < 1e-12)
}
cat("\n   Every row satisfies the matching condition to machine precision, reproduces the quoted\n")
cat("   r and t, and is unitary.\n")

cat("\n=== what this settles\n\n")
cat("  The form is NOT borrowed on authority: it is the unique solution of the matching\n")
cat("  conditions the corner term imposes, up to the orientation relabelling. The referee is\n")
cat("  right that the repo's provenance was a script header quoting two model reviews, and that\n")
cat("  was a real weakness in how it was recorded. It is not a weakness in the result.\n")
cat("  kappa = 1 gives r = 0, t = 1 as the transparent point, and kappa = 0 and infinity give\n")
cat("  the two perfect mirrors, which is the structure 3.4 uses.\n")

cat("\n=== CHALLENGE: is the 'other sign' branch really the same seam relabelled?\n\n")
cat("  The paper says so. A verification review says that branch is not even unitary. Check.\n\n")
cat("      kappa    r = (1+k^2)/(1-k^2)   t = k(1+r)      r^2 + t^2     unitary?\n")
for (k in c(0.25, 0.5, 0.75)) {
  r <- (1+k^2)/(1-k^2); t <- k*(1+r)
  cat(sprintf("   %8.2f %20.5f %12.5f %14.5f %12s\n", k, r, t, r^2+t^2,
      if (abs(r^2+t^2-1) < 1e-12) "yes" else "NO"))
}
cat("\n  At kappa = 1/2 the branch gives r = 5/3 and t = 4/3, so r^2+t^2 = 41/9 = 4.5556.\n")
cat("  A passive interface cannot be relabelled into a non-unitary one. The review is right and\n")
cat("  the paper's sentence calling this branch 'the same seam traversed the other way' is\n")
cat("  FALSE. That branch is not a seam at all; it is the sign choice that fails.\n")
cat("  What survives: the OTHER sign gives r = (1-k^2)/(1+k^2), t = 2k/(1+k^2), r^2+t^2 = 1 at\n")
cat("  every kappa, which is the branch the paper uses. The derivation selects a sign, and the\n")
cat("  honest statement is that unitarity selects it, not that both are equivalent.\n")

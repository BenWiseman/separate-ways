# A.13 argues that a constant reflectivity puts an image term in the two-point function which is
# singular between mirror-placed points on OPPOSITE sides of the seam, and that no Hadamard state
# is singular there, so r = 0. The verification pass says the cross-side image term cancels for a
# unitary two-port defect and the argument is therefore empty. Check it.
#
# Scattering states on a line with a point defect, S = [[r,t],[t,-r]], r^2 + t^2 = 1.
#   left-incident :  x<0 : e^{ikx} + r e^{-ikx}        x>0 : t e^{ikx}
#   right-incident:  x>0 : e^{-ikx} + r' e^{ikx}       x<0 : t' e^{-ikx}     with r' = -r, t' = t
# The mode sum for x<0<y collects a coefficient on e^{-ik(x+y)}, which is the mirror-coincidence
# structure A.13 needs. Compute that coefficient.
coef_cross <- function(r) {
  t <- sqrt(1-r^2); rp <- -r; tp <- t
  # left channel  : (e^{ikx} + r e^{-ikx}) * conj(t e^{iky})  -> r*Conj(t) on e^{-ik(x+y)}
  # right channel : (t' e^{-ikx}) * conj(e^{-iky} + r' e^{iky}) -> t'*Conj(r') on e^{-ik(x+y)}
  r*Conj(t) + tp*Conj(rp)
}
# and the SAME-side coefficient, for contrast: both channels contribute on x,y < 0
coef_same <- function(r) {
  t <- sqrt(1-r^2); tp <- t
  # left channel gives r on e^{-ik(x+y)}; right channel contributes |t'|^2 only to e^{-ik(x-y)}
  r
}
cat("      r        t        cross-side image coeff     same-side image coeff\n")
for (r in c(0.1, 0.25, 0.5, 0.75, 0.9, 0.99)) {
  cat(sprintf("  %6.3f  %7.4f  %22.3e  %20.4f\n",
      r, sqrt(1-r^2), Mod(coef_cross(r)), coef_same(r)))
}
cat("\n  The cross-side coefficient is r*t + t*(-r) = 0 identically: unitarity of the two-port\n")
cat("  defect cancels it. The mirror-coincidence singularity A.13 forbids never appears ACROSS\n")
cat("  the seam, so forbidding it constrains nothing and cannot force r = 0.\n")
cat("\n  The same-side image term is r and does NOT cancel, which is the structure that genuinely\n")
cat("  exists. It is singular at x = -y with both points on the same side. Whether THAT violates\n")
cat("  Hadamard regularity is a separate question this file does not settle, and it is the\n")
cat("  question A.13 would have to answer to recover a derivation.\n")
cat("\n  Consistent with Mintchev-Sorba: reflection terms have same-side support, transmission\n")
cat("  terms have opposite-side support. A.13 placed a reflection term on opposite sides.\n")

cat("\n  WHY THE ARGUMENT FAILS COMPLETELY, rather than just on one side.\n")
cat("  Write the mirror of y as ybar = -y. The image term is singular when x = ybar, i.e. x + y = 0.\n")
for (case in c("same side (x<0, y<0)", "opposite sides (x<0, y>0)")) {
  if (grepl("same", case)) {
    xs <- c(-0.1,-1,-5); ys <- c(-0.1,-2,-7)
    cat(sprintf("   %s : image coefficient = r (survives)\n", case))
    cat("      but x + y is a sum of two negatives, so x + y = 0 is unreachable:\n")
    for (i in seq_along(xs)) cat(sprintf("        x=%5.1f y=%5.1f  x+y=%6.1f\n", xs[i], ys[i], xs[i]+ys[i]))
  } else {
    cat(sprintf("   %s : x + y = 0 IS reachable (e.g. x=-3, y=3)\n", case))
    cat("      but the image coefficient there is identically 0, as the table above shows.\n")
  }
}
cat("\n  So where the image term exists, its singularity cannot be reached; where the singularity\n")
cat("  could be reached, the term is absent. The inference has no case left to bite on, and\n")
cat("  Hadamard regularity does not force r = 0. Transparency has to be ADOPTED, not derived.\n")

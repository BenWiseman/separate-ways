# Ben: is Hawking radiation needed at all, if matter transits the singularity into a past
# white hole? The fold's own machinery answers this and the answer is not the hoped one.
#
# Last night: alpha^2 = 1 reads W(t - i beta) = W(t), which HOLDS ONLY at beta = 2 pi/H
# (2 pi / kappa for a black hole). That is the KMS condition at the Hawking temperature.
# So the fold's involution property does not merely permit the Hawking temperature, it
# REQUIRES it. Quantify: how far from thermal can a state be and still satisfy it?
#
# KMS in frequency space is detailed balance:  Wt(-w) = e^{-beta w} Wt(w).
# Build a state whose occupation departs from the Bose factor by a parameter and measure
# the violation.

beta <- 1
nB <- function(w) 1/(exp(beta*w)-1)                  # thermal occupation
# departure: n(w) = nB(w) * (1 + d * f(w)), f bounded, d the departure parameter
Wt <- function(w, d, f) { n <- nB(abs(w))*(1 + d*f(abs(w)))
  ifelse(w > 0, 1 + n, n) }                          # emission vs absorption weights

cat("=== 1. exact thermality gives detailed balance; anything else does not\n\n")
cat("   R(w) = Wt(-w) / (e^{-beta w} Wt(w)).  R = 1 is the KMS condition.\n\n")
f <- function(w) w/(1+w)                             # a smooth bounded departure
cat("        d         R(w=0.5)     R(w=1.5)     R(w=3.0)    max |R-1|\n")
for (d in c(0, 1e-4, 1e-2, 0.1, 0.5)) {
  ws <- c(0.5, 1.5, 3.0)
  R <- sapply(ws, function(w) Wt(-w,d,f)/(exp(-beta*w)*Wt(w,d,f)))
  Rall <- sapply(seq(0.05,6,length.out=400), function(w) Wt(-w,d,f)/(exp(-beta*w)*Wt(w,d,f)))
  cat(sprintf("   %8.0e %12.6f %12.6f %12.6f %12.3e\n", d, R[1],R[2],R[3], max(abs(Rall-1))))
}
cat("\n  Only the exactly thermal occupation satisfies it. The violation is first order\n")
cat("  in the departure, so there is no window in which a nearly-thermal state passes.\n")

cat("\n=== 2. what that means for the three things Ben's question bundles together\n\n")
cat("  (a) The Hawking TEMPERATURE. Forced. alpha^2 = 1 is exactly the statement that\n")
cat("      the state is KMS at 2 pi/kappa, and last night's sweep showed it fails by\n")
cat("      O(25) at every other period. The fold cannot drop the temperature: its own\n")
cat("      involution property is the thing that fixes it. Dropping it breaks the fold.\n\n")
cat("  (b) The Hawking FLUX. Not forced, and not forbidden, because it is outside the\n")
cat("      domain. A state in equilibrium at T_H (Hartle-Hawking) is KMS and has no NET\n")
cat("      flux. An evaporating hole is in the Unruh state, which is NOT KMS in Killing\n")
cat("      time, so alpha^2 = 1 does not hold for it and none of A.10 to A.14 applies.\n")
cat("      That is the same domain boundary A.14 already records for collapse holes,\n")
cat("      arrived at from a different direction, which is a consistency check on both.\n\n")
cat("  (c) Transit into a past white hole as a REPLACEMENT for evaporation. It cannot\n")
cat("      be one. The image of infalling matter sits in P, which is in the PAST of both\n")
cat("      exteriors, so nothing is returned to our future by that route. Evaporation is\n")
cat("      a statement about our future. The two are not alternatives; they are answers\n")
cat("      to different questions.\n")

cat("\n=== 3. flatly\n\n")
cat("  The route 'drop Hawking radiation because matter transits the singularity' is\n")
cat("  closed, and closed by the fold's own machinery rather than by deference to the\n")
cat("  literature. The temperature is forced by alpha^2 = 1. The transit does not point\n")
cat("  at our future and so cannot do evaporation's job.\n")

cat("\n=== 4. the route, same turn\n\n")
cat("  What the calculation does open is a sharper domain statement than A.14 has. The\n")
cat("  paper's horizon machinery describes the EQUILIBRIUM state, and the object it\n")
cat("  describes is an eternal hole in the Hartle-Hawking state, not an astrophysical\n")
cat("  hole that formed and will evaporate. Two consequences worth chasing:\n")
cat("   - The cosmological horizon IS in equilibrium and is genuinely KMS, so it is the\n")
cat("     one horizon where all of A.10 to A.14 applies without apology. That agrees\n")
cat("     with A.14's own domain paragraph, reached independently.\n")
cat("   - The size of the KMS violation for a slowly evaporating hole is computable and\n")
cat("     small: the Unruh state departs from Hartle-Hawking by O(1/(kappa t_evap)). If\n")
cat("     the fold's predictions degrade smoothly in that parameter rather than failing\n")
cat("     outright, the machinery extends to astrophysical holes as a controlled\n")
cat("     approximation. That is the calculation that would widen the domain, and it is\n")
cat("     not done here.\n")

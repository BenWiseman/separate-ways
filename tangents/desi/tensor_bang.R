# The consequence 4.2 states and never collects.
#
# Tensor perturbations on FRW obey, for v = a h,   v'' + (k^2 - a''/a) v = 0.
# That is the minimally coupled massless equation. Section 4.2 already records that at a
# radiation bang a'' = 0. So the source term for gravitons vanishes identically and the
# bang produces NO tensor background. The paper uses this fact to argue that massless
# fields cannot distinguish the two orientations, and then drops it.

cat("=== 1. a'' at a bang, as a function of the equation of state\n\n")
cat("   a  proportional to  eta^q  with  q = 2/(1+3w),  so  a''/a = q(q-1)/eta^2.\n\n")
cat("        w        q = 2/(1+3w)     q(q-1)      graviton source\n")
for (w in c(1, 1/3, 0.30, 0.25, 0, -1/3)) {
  q <- 2/(1+3*w)
  cat(sprintf("   %8.4f %14.4f %12.5f %18s\n", w, q, q*(q-1),
      if (abs(q*(q-1)) < 1e-12) "EXACTLY ZERO" else "nonzero"))
}
cat("\n   Radiation is the one equation of state that switches the source off exactly.\n")
cat("   It is also the one BFT's bang has. That is not a coincidence the paper chose;\n")
cat("   A.5 shows radiation is forced as the saturating case at the bang.\n")

cat("\n=== 2. how fast the source returns as w departs from a third\n\n")
cat("   Writing w = 1/3 - e:   q(q-1) = 3e/2 + O(e^2), so the source is FIRST ORDER in\n")
cat("   the departure and the produced energy, going as amplitude squared, is second.\n\n")
cat("      departure e     q(q-1)        suppression vs matter (e = 1/3)\n")
for (e in c(1e-4, 1e-3, 1e-2, 1e-1, 1/3)) {
  w <- 1/3 - e; q <- 2/(1+3*w); s <- (q*(q-1))^2 / (2^2)
  cat(sprintf("   %12.0e %12.5f %30.3e\n", e, q*(q-1), s))
}
cat("\n   A one per cent departure from pure radiation suppresses the tensor energy by\n")
cat("   about four orders of magnitude relative to a matter era. At the bang every\n")
cat("   species is relativistic and the departure is far below that.\n")

cat("\n=== 3. the Hadamard tail for a massless field: no floor to hide behind\n\n")
cat("   Section 5.1's ultraviolet seed is c_+(0)/c_-(0) ~ i*gamma/(4 p^2) with gamma = m a_1.\n")
cat("   For the graviton m = 0, so gamma = 0 and the seed VANISHES. The allowed power-law\n")
cat("   part of the occupation is therefore identically zero, and the whole graviton\n")
cat("   occupation must fall faster than any power of k.\n")
cat("   For the massive field the seed is nonzero, which is exactly why A.1's p^-4\n")
cat("   bang-adiabatic tail is a near miss rather than nonsense: it has the right power\n")
cat("   but the wrong coefficient, and its energy integral diverges logarithmically.\n")
cat("   That divergence is the m=1 moment of the general condition. The paper already\n")
cat("   uses it to exclude a rival state without naming it as one case of a criterion.\n")

cat("\n=== 4. the falsifier, stated so it can fail\n\n")
cat("   BFT's bang makes no primordial tensor background. Inflation makes one at every\n")
cat("   frequency. So ANY stochastic gravitational-wave background established as\n")
cat("   primordial rather than astrophysical, at ANY frequency, falsifies the bang this\n")
cat("   paper's section 5 rests on. No amplitude needs to be predicted for the test to\n")
cat("   bite, which is unusual and makes it cheap to state and hard to wriggle out of.\n")

cat("\n=== 5. and this is where DESI enters, which is the question asked\n\n")
cat("   The live case is the nanohertz background pulsar timing arrays have detected.\n")
cat("   Two readings compete: supermassive black-hole binaries, or something primordial.\n")
cat("   Under the fold the second reading is excluded, so the fold PREDICTS the nanohertz\n")
cat("   background is entirely astrophysical.\n")
cat("   DESI measures the galaxy population and merger history that sets the astrophysical\n")
cat("   amplitude. If DESI's population cannot supply the observed amplitude, the shortfall\n")
cat("   is primordial and the bang is falsified. That is a spectroscopy-times-gravitational\n")
cat("   -wave test of this paper, on data that already exists.\n")
cat("   It also couples to 5.4: little red dots, if they are black holes, RAISE the\n")
cat("   astrophysical budget and so help the fold rather than embarrass it.\n")

# The shape of the fold term across the contact surface. Ben asked whether the coefficient could
# be a relu or a sigmoid rather than a number. It is closer to the first than to either.
#
# swept_length.R fixes the amplitude AT the boundary. This asks what the term does either side.

M <- 1
slope <- 3 * pi + 8            # sigma -> (3 pi + 8) M (r - M), from the shooting solution
A <- 3.9004                    # Delta^{1/2} -> A M s^{-1/2}, from swept_length.R

cat("=== 1. the three regions, from the causal structure alone ===\n")
cat("   contact_world_function.py established sigma(x, Theta x) = (3 pi + 8) M (r - M) near the\n")
cat("   boundary: positive outside, zero at r = M, negative inside. Those are three different\n")
cat("   physical situations and the image term does something different in each.\n\n")
cat("      r/M      sigma            relation to the image      image term\n")
for (r in c(1.06, 1.02, 1.0, 0.98, 0.94)) {
  s <- slope * M * (r - M)
  rel <- if (s > 1e-12) "spacelike" else if (s < -1e-12) "timelike" else "NULL"
  beh <- if (s > 1e-12) "finite and real" else if (s < -1e-12) "complex" else "divergent"
  cat(sprintf("   %7.2f  %+14.6f   %-24s  %s\n", r, s, rel, beh))
}

cat("\n=== 2. outside: not small, exactly zero ===\n")
cat("   Beyond the contact boundary the pair is spacelike and the image term is finite, and\n")
cat("   beyond the HORIZON the silence theorem puts the cross-sheet commutator at zero\n")
cat("   identically. So the term does not taper off, it stops. That is the relu-like part of\n")
cat("   the shape, and it is a theorem rather than a suppression.\n")

cat("\n=== 3. at the boundary: a one-sided divergence with a known power and now a number ===\n")
cat("   G_img ~ A delta^{-3/2} with delta the distance off r = M, and A fixed by the swept\n")
cat("   length. In terms of sigma the divergence is sigma^{-3/2}, so\n")
for (d in c(1e-1, 1e-2, 1e-3)) {
  cat(sprintf("      delta = %.0e M:  sigma = %.4e,  delta^{-3/2} = %10.1f,  A delta^{-3/2} = %12.1f\n",
              d, slope * d, d^(-1.5), A * d^(-1.5)))
}

cat("\n=== 4. inside: the term goes complex, which is the part worth having ===\n")
cat("   For sigma < 0 the pair is timelike separated and sigma^{-3/2} is not real. Continuing\n")
cat("   below the cut, sigma^{-3/2} -> |sigma|^{-3/2} e^{-3 pi i / 2} = i |sigma|^{-3/2}, so the\n")
cat("   image term acquires an IMAGINARY part inside the contact surface and the real part\n")
cat("   vanishes there:\n\n")
cat("      r/M      |sigma|        Re part        Im part\n")
for (r in c(0.99, 0.96, 0.90, 0.80)) {
  s <- slope * M * (r - M); z <- A * as.complex(s)^(-1.5)
  cat(sprintf("   %7.2f  %12.6f  %13.6f  %13.6f\n", r, abs(s), Re(z), Im(z)))
}
cat("\n   An imaginary part in the two-point function is a decay rate. So the fold's contact\n")
cat("   region is not merely a place where a stress diverges: it is a region where the\n")
cat("   effective action has an imaginary part, which means pair creation.\n")
cat("   That is a statement the construction had not made, and it comes free with the sign of\n")
cat("   sigma rather than from any new assumption.\n")

cat("\n=== 5. so the shape, in one line each ===\n")
cat("      outside the horizon      exactly zero, by the silence theorem\n")
cat("      horizon to r = M         finite, real, falling as delta^{-3/2} towards the boundary\n")
cat("      at r = M                 divergent, coefficient 3.90 M\n")
cat("      inside r = M             imaginary, which is a decay rate\n")
cat("   Relu-like at the outer edge because the cut-off is a theorem; a pole at the boundary;\n")
cat("   and a change of character rather than of size on the inside. Nothing about it is a\n")
cat("   sigmoid, and nothing about it is a constant.\n")

cat("\n=== 6. the plant ===\n")
cat("   The continuation must give a real answer outside and a complex one inside, or it is not\n")
cat("   tracking the causal structure at all.\n")
for (s in c(0.5, -0.5)) {
  z <- as.complex(s)^(-1.5)
  cat(sprintf("      sigma = %+0.1f:  sigma^{-3/2} = %s   %s\n", s, format(z, digits = 6),
              ifelse(abs(Im(z)) < 1e-12, "real, as it must be outside", "complex, as it must be inside")))
}
stopifnot(abs(Im(as.complex(0.5)^(-1.5))) < 1e-12, abs(Im(as.complex(-0.5)^(-1.5))) > 0.1)

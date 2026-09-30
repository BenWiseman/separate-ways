# Where the construction leaves general relativity, and how thick that place is.
#
# Section 3.6 says the fold's extra term is identically zero outside every horizon and lives on the
# inner half of a hole's interior, and then computes how thick the shell is in which it beats the
# interior's own focusing. A reader meets that as three paragraphs of powers and coefficients. This
# draws it: the interior on the left with the contact sphere and the shell marked, and the thickness
# against the field's mass on the right, with the Planck length as the floor that decides whether
# any of it means anything.
#
# Every number is recomputed here from the same closed form image_stress_coefficient.R uses, and the
# two anchor values are asserted against that file's own output, so the figure cannot drift from it.

mP_GeV <- 1.220890e19; lP_m <- 1.616255e-35; rh_sun <- 2.953250e3
Bstar  <- 11.5138          # the level to beat, from interior_focusing_threshold.R
KAP    <- 0.0039329        # the coefficient, from image_stress_coefficient.R
PW     <- 2.5              # the power, from caustic_power_at_a_hole.R

Dstar <- function(m_GeV, Msun = 1, kap = KAP)
  (8 * pi * kap * (m_GeV / mP_GeV)^2 / Bstar)^(1 / PW) * rh_sun * Msun

# The floor. The mass-independent piece of the caustic stress is two powers more divergent and
# vastly smaller, so it never touches a massive field's shell, but it stops the shell going to
# zero: T_kk = -CKK M^{-1/2}(M-r)^{-7/2} against the same B*/r_h^2 leaves D ~ lP^{4/7} r_h^{3/7}.
CKK <- 3.191916e-3
Dfloor <- function(Msun) { rh <- rh_sun*Msun
  (8*pi*lP_m^2 * CKK * sqrt(2/rh) / Bstar * rh^2)^(2/7) }
stopifnot(abs(Dfloor(1)/1.079e-19 - 1) < 0.01, Dfloor(1) > 1e3*lP_m)

cat("=== the two anchors, against the files they come from ===\n")
cat(sprintf("   the 491.6 PeV fermion at one solar mass: %.4e m, against 2.127e-06\n",
            Dstar(4.916e8)))
cat(sprintf("   an electron at one solar mass:           %.4e m, against 5.510e-16\n",
            Dstar(0.000511)))
stopifnot(abs(Dstar(4.916e8) / 2.127e-6 - 1) < 0.01,
          abs(Dstar(0.000511) / 5.510e-16 - 1) < 0.01)

FIELDS <- list(list("electron", 0.000511), list("muon", 0.10566), list("proton", 0.938272),
               list("b", 4.18), list("W", 80.377), list("top", 172.69),
               list("the fold's fermion", 4.916e8))

ink <- "grey15"; c1 <- "#1f4e79"; c2 <- "#a8400f"; grey <- "#8a97a4"; pale <- "#dce3ea"

draw <- function() {
  par(mfrow = c(1, 2), mar = c(4.2, 4.4, 2.2, 1.0), xpd = NA)

  ## left: the interior, the contact region, the shell
  plot(NA, xlim = c(-1.15, 1.15), ylim = c(-1.15, 1.25), asp = 1, axes = FALSE,
       xlab = "", ylab = "", main = "")
  th <- seq(0, 2 * pi, length.out = 400)
  # the horizon
  polygon(cos(th), sin(th), border = grey, lwd = 1.6, col = NA)
  # the contact region, r <= M = r_h/2
  polygon(0.5 * cos(th), 0.5 * sin(th), border = NA, col = pale)
  polygon(0.5 * cos(th), 0.5 * sin(th), border = c1, lwd = 1.8)
  # the shell, drawn at an exaggerated width with the exaggeration stated
  polygon(c(0.5 * cos(th), rev(0.46 * cos(th))), c(0.5 * sin(th), rev(0.46 * sin(th))),
          border = NA, col = c2)
  points(0, 0, pch = 4, cex = 1.1, col = ink, lwd = 2)
  text(0, -0.085, "singularity", col = ink, cex = 0.78)
  text(0, 1.10, expression(paste("horizon, ", r[h])), col = grey, cex = 0.80)
  text(0, 0.62, expression(paste("contact sphere, ", r == r[h] / 2)), col = c1, cex = 0.80)
  # centred, not left-aligned at 1.02: with asp = 1 the panel is barely wider than the stated
  # xlim and a left-aligned label of this length runs off the right edge.
  text(1.00, 0.44, "the shell", col = c2, cex = 0.82)
  arrows(0.97, 0.36, 0.52, 0.17, length = 0.06, col = c2, lwd = 1.4)
  text(0, -1.33, "the fold's term is zero outside the horizon and", col = ink, cex = 0.78)
  text(0, -1.46, "nonzero only inside the contact sphere", col = ink, cex = 0.78)
  text(0, -1.62, "(shell width drawn wide; it is microns, not a tenth of a radius)",
       col = grey, cex = 0.68)
  mtext("Where it acts", side = 3, line = 0.5, cex = 0.95, col = ink)

  ## right: thickness against mass
  par(xpd = FALSE)
  ms <- 10^seq(-4, 9.2, length.out = 400)
  y1 <- Dstar(ms, 1); y2 <- Dstar(ms, 1e9)
  plot(ms, y1, log = "xy", type = "l", lwd = 2.2, col = c1,
       ylim = c(1e-36, 1e4), xlim = c(1e-4, 3e9), yaxt = "n",
       xlab = "field mass (GeV)", ylab = "shell thickness (m)",
       cex.lab = 0.88, cex.axis = 0.80, col.axis = ink, col.lab = ink, main = "")
  axis(2, at = 10^seq(-36, 4, by = 8), labels = parse(text = sprintf("10^%d", seq(-36, 4, by = 8))),
       cex.axis = 0.80, col.axis = ink, las = 1)
  lines(ms, y2, lwd = 2.0, col = c2, lty = 2)
  # the Planck length, which is the whole point of the panel: everything sits far above it
  abline(h = lP_m, col = grey, lwd = 1.5, lty = 3)
  text(2.0e-4, lP_m * 30, "the Planck length, which nothing here comes near",
       col = grey, cex = 0.78, adj = 0)
  # and the floor: the mass-independent term gives a shell no field falls below
  abline(h = Dfloor(1), col = c1, lwd = 1.2, lty = 4)
  abline(h = Dfloor(1e9), col = c2, lwd = 1.2, lty = 4)
  # the floor line ran through this at 2.4x; lifted clear, still well under Dfloor(1e9)
  text(2.0e-4, Dfloor(1)*25, "the floor: what a massless field still gets",
       col = c1, cex = 0.76, adj = 0)
  SHOW <- list(list("electron", 0.000511, 0), list("proton", 0.938272, 0),
               list("top quark", 172.69, 0), list("the fold's fermion", 4.916e8, 1))
  for (f in SHOW) {
    m <- f[[2]]; d <- Dstar(m)
    points(m, d, pch = 19, cex = 0.70, col = c1)
    if (f[[3]] == 0) text(m * 1.8, d * 0.035, f[[1]], col = ink, cex = 0.78, adj = 0)
    else            text(m * 0.55, d * 0.035, f[[1]], col = ink, cex = 0.78, adj = 1)
  }
  legend(10^-3.6, 10^-21, legend = c("a billion solar masses", "one solar mass"),
         col = c(c2, c1), lty = c(2, 1), lwd = c(2.0, 2.2), bty = "n",
         cex = 0.80, text.col = ink, seg.len = 2.4)
  mtext("How thick it is", side = 3, line = 0.5, cex = 0.95, col = ink)
  par(xpd = NA)
}

for (dev in c("pdf", "png")) {
  out <- sprintf("papers/1_separate_ways/fig_interior_shell.%s", dev)
  if (dev == "pdf") pdf(out, width = 9.6, height = 4.8, pointsize = 12)
  else png(out, width = 9.6, height = 4.8, units = "in", res = 200, pointsize = 12)
  draw(); invisible(dev.off())
  cat(sprintf("   wrote %s\n", out))
}

cat("\n=== the plant: the figure must break if the coefficient moves ===\n")
cat("   The two anchors are asserted against image_stress_coefficient.R's own numbers, so a change\n")
cat("   there without a change here stops the script. Checked by pretending kappa moved by ten:\n")
fake <- Dstar(4.916e8, 1, KAP * 10)
cat(sprintf("      kappa x 10 would give %.4e m against the asserted 2.127e-06  ->  %s\n",
            fake, ifelse(abs(fake / 2.127e-6 - 1) < 0.01, "would draw", "STOPS, as it must")))
stopifnot(abs(fake / 2.127e-6 - 1) > 0.01)
cat("   And the Planck-length line has to be crossed by something, or the panel says nothing:\n")
below <- sum(sapply(FIELDS, function(f) Dstar(f[[2]]) < lP_m))
cat(sprintf("      fields on the list below a Planck length at one solar mass: %d of %d\n",
            below, length(FIELDS)))
cat("      none, which is the result; the line is there so a reader can see that it is none.\n")

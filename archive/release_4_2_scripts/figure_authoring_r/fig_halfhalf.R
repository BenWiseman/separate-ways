# FIGURE: "What the fold does to one mode at the bang."
#
# Section 2.2's result -- Theta-invariance forces |psi_1(0)| = |psi_2(0)| = 1/sqrt(2), so every
# mode is exactly half in each component at the bang, and what survives is one phase per mode
# -- is reached through sigma_x, a starred wavefunction and "invariance as a ray". None of
# those is a picture. This draws the two statements instead.
#
# The curves are an honest RK4 integration of the mode equation, not a sketch, and the script
# checks the Theta relation it is meant to illustrate before drawing anything.
#
# Base R only. Authored at 6.5 in, the width the manuscript places it at.

INK <- "#1b1b1b"; C1 <- "#0b6e9e"; C2 <- "#c2461f"; GREY <- "#7e7e7e"
FAM <- c("#0b6e9e", "#8a6d1f", "#5f7a3a")

GAMMA <- 1.0; P <- 0.60

# i dpsi/deta = H psi,  H = [[gamma*eta, p], [p, -gamma*eta]]
Hmat <- function(e) matrix(c(GAMMA * e, P, P, -GAMMA * e), 2, 2, byrow = TRUE)
rhs  <- function(e, y) as.vector(-1i * (Hmat(e) %*% y))

rk4 <- function(psi0, eta1, n = 4000) {
  h <- eta1 / n
  etas <- numeric(n + 1); out <- matrix(0 + 0i, n + 1, 2)
  psi <- psi0; eta <- 0
  etas[1] <- 0; out[1, ] <- psi
  for (i in 1:n) {
    k1 <- rhs(eta,       psi)
    k2 <- rhs(eta + h/2, psi + (h/2) * k1)
    k3 <- rhs(eta + h/2, psi + (h/2) * k2)
    k4 <- rhs(eta + h,   psi + h * k3)
    psi <- psi + (h/6) * (k1 + 2*k2 + 2*k3 + k4)
    eta <- eta + h
    etas[i + 1] <- eta; out[i + 1, ] <- psi
  }
  list(eta = etas, psi = out)
}

# a Theta-invariant solution: psi(0) = (1, e^{i mu}) / sqrt(2)
solve_mode <- function(mu, T = 4) {
  psi0 <- c(1 + 0i, exp(1i * mu)) / sqrt(2)
  fw <- rk4(psi0,  T); bw <- rk4(psi0, -T)
  ord <- rev(seq_along(bw$eta))
  list(eta = c(bw$eta[ord], fw$eta[-1]),
       psi = rbind(bw$psi[ord, , drop = FALSE], fw$psi[-1, , drop = FALSE]))
}

MU_SET <- c(0, 2*pi/3, 4*pi/3)
SOL <- lapply(MU_SET, solve_mode)

# ---- checks. A figure that illustrates a theorem should verify it first.
for (k in seq_along(SOL)) {
  s <- SOL[[k]]; mu <- MU_SET[k]
  nrm <- rowSums(Mod(s$psi)^2)
  stopifnot(max(abs(nrm - 1)) < 1e-9)                       # unitary evolution
  i0 <- which.min(abs(s$eta))
  stopifnot(abs(Mod(s$psi[i0, 1])^2 - 0.5) < 1e-12,         # half in each component
            abs(Mod(s$psi[i0, 2])^2 - 0.5) < 1e-12)
  # Theta: sigma_x conj(psi(-eta)) == e^{-i mu} psi(eta)
  n <- length(s$eta); flip <- n:1
  lhs <- cbind(Conj(s$psi[flip, 2]), Conj(s$psi[flip, 1]))
  stopifnot(max(Mod(lhs - exp(-1i * mu) * s$psi)) < 1e-8)
}
cat("checks passed: norm, half-and-half at the bang, and the Theta relation\n")

draw <- function() {
  par(mfrow = c(1, 2), mar = c(3.7, 4.1, 2.3, 1.0), family = "sans")
  YL <- c(-0.08, 1.48)

  # guides drawn as segments, NOT abline: with xpd on, abline runs the whole device width
  guides <- function() {
    segments(0, -0.02, 0, 1.06, col = GREY, lwd = 1.2, lty = 2)
    segments(-4, 0.5, 4, 0.5, col = GREY, lwd = 1.0, lty = 3)
  }
  frame_axes <- function(ylab = NULL) {
    axis(1, at = seq(-4, 4, 2), cex.axis = 0.92, col = GREY, col.axis = INK,
         tck = -0.02, mgp = c(2, 0.6, 0))
    axis(2, at = c(0, 0.5, 1), labels = c("0", "1/2", "1"), las = 1, cex.axis = 0.92,
         col = GREY, col.axis = INK, tck = -0.02, mgp = c(2, 0.7, 0))
    mtext("time, through the bang", side = 1, line = 2.3,
          cex = 0.90, col = INK)
    if (!is.null(ylab)) mtext(ylab, side = 2, line = 2.5, cex = 0.90, col = INK)
  }

  # ---------------- (a) the constraint ----------------
  s1 <- SOL[[1]]
  plot(NA, xlim = c(-4, 4), ylim = YL, axes = FALSE, xlab = "", ylab = "")
  mtext("(a)  exactly half, at the bang", side = 3, line = 0.6, adj = 0,
        cex = 0.98, font = 2, col = INK)
  guides()
  lines(s1$eta, Mod(s1$psi[, 1])^2, col = C1, lwd = 2.7)
  lines(s1$eta, Mod(s1$psi[, 2])^2, col = C2, lwd = 2.7)
  points(0, 0.5, pch = 19, cex = 1.3, col = INK)
  frame_axes("share of the mode in each component")

  # legend in the clear band above the curves, so no label sits on a curve
  segments(-3.9, 1.40, -3.3, 1.40, col = C1, lwd = 2.7)
  text(-3.15, 1.40, "one component", col = INK, cex = 0.90, adj = 0)
  segments(-3.9, 1.25, -3.3, 1.25, col = C2, lwd = 2.7)
  text(-3.15, 1.25, "the other", col = INK, cex = 0.90, adj = 0)
  points(-3.6, 1.10, pch = 19, cex = 1.2, col = INK)
  text(-3.15, 1.10, "the fold pins this crossing", col = INK, cex = 0.90, adj = 0)

  # ---------------- (b) the surviving freedom ----------------
  plot(NA, xlim = c(-4, 4), ylim = YL, axes = FALSE, xlab = "", ylab = "")
  mtext("(b)  one phase survives", side = 3, line = 0.6, adj = 0,
        cex = 0.98, font = 2, col = INK)
  guides()
  for (k in seq_along(SOL)) lines(SOL[[k]]$eta, Mod(SOL[[k]]$psi[, 1])^2,
                                  col = FAM[k], lwd = 2.4)
  points(0, 0.5, pch = 19, cex = 1.3, col = INK)
  frame_axes()

  text(-3.9, 1.40, "three members of the family:", col = INK, cex = 0.90, adj = 0)
  text(-3.9, 1.25, "pinned there, free elsewhere", col = INK, cex = 0.90, adj = 0)
}

dir.create("pub/paper2/figs", showWarnings = FALSE, recursive = TRUE)
pdf("pub/paper2/figs/fig_halfhalf.pdf", width = 6.5, height = 3.5); draw(); dev.off()
png("pub/paper2/figs/fig_halfhalf.png", width = 1950, height = 1050, res = 300); draw(); dev.off()
cat("wrote fig_halfhalf.{pdf,png}\n")

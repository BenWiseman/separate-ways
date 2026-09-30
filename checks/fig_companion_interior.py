"""The contact region drawn as concentric bands through a Schwarzschild interior.

Colour carries a quantity and not a mood. Every band is shaded by

    dphi_max(r) = 2 pi - 4 arcsin sqrt(r / 2M),

the angle a causal curve starting at radius r can still sweep on its two interior legs
before it reaches the singularity. The antipodal map puts a point half a turn from its
own image, so the sheets touch wherever that angle reaches pi, which happens at r = M and
holds everywhere inside. Bands outside r = M are held in a flat sepia and bands inside it
carry the spectrum, so the colour starts exactly where the contact region does.

Checks at the bottom: the three anchor values, the crossing radius solved independently
by bisection rather than read off the closed form, and a deliberately wrong exponent that
must fail.

Writes fig_companion_interior.pdf and .png next to this file.
"""
import colorsys
import math
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Wedge

OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                   "pub", "paper2")
M = 1.0
RH = 2.0 * M                      # horizon radius in these units
NBANDS = 420


def dphi_max(r):
    """Angle available to a causal curve on its two interior legs, from radius r."""
    return 2.0 * math.pi - 4.0 * math.asin(math.sqrt(r / RH))


def band_colour(r):
    """Sepia outside the contact radius, spectrum inside it."""
    if r > M:
        u = (RH - r) / (RH - M)                     # 0 at the horizon, 1 at r = M
        return colorsys.hls_to_rgb(36 / 360.0, 0.90 - 0.13 * u, 0.20)
    t = (M - r) / M                                 # 0 at r = M, 1 at the singularity
    return colorsys.hls_to_rgb((6 + 288 * t) / 360.0, 0.80 - 0.05 * t, 0.62)


# Sixteen bands at equal steps of the angle, so every band is one sixteenth of a turn
# and r = M falls exactly on the eighth contour. Eight sepia below half a turn, eight
# spectral above it. r(dphi) inverts dphi_max.
NB = 16
EDGES = [RH * math.sin((2 * math.pi - k * 2 * math.pi / NB) / 4.0) ** 2
         for k in range(NB + 1)]                     # EDGES[0] = 2M ... EDGES[NB] = 0


def step_colour(k):
    """Colour of the band between EDGES[k] and EDGES[k+1]; k = 0 is just inside the horizon."""
    if k < NB // 2:
        u = k / (NB // 2 - 1.0)
        return colorsys.hls_to_rgb(36 / 360.0, 0.915 - 0.135 * u, 0.22)
    t = (k - NB // 2) / (NB // 2 - 1.0)
    return colorsys.hls_to_rgb((4 + 292 * t) / 360.0, 0.815 - 0.055 * t, 0.66)


fig = plt.figure(figsize=(7.4, 3.6))
fig.patch.set_facecolor("white")
gs = fig.add_gridspec(1, 2, width_ratios=[1.18, 1.0], wspace=0.26,
                      left=0.025, right=0.975, top=0.795, bottom=0.155)

# ----------------------------------------------------------------- left: the map
ax = fig.add_subplot(gs[0])
for k in range(NB):
    ax.add_patch(Wedge((0, 0), EDGES[k], 0, 180, facecolor=step_colour(k),
                       edgecolor="white", linewidth=0.55, zorder=k + 1))

th = np.linspace(0, math.pi, 400)
ax.plot(RH * np.cos(th), RH * np.sin(th), color="#2b2b2b", lw=1.6, zorder=40)
ax.plot([-RH, RH], [0, 0], color="#2b2b2b", lw=1.6, zorder=40)
ax.plot(M * np.cos(th), M * np.sin(th), color="#3a3a3a", lw=1.5,
        ls=(0, (4.5, 2.6)), zorder=41)

ax.annotate("horizon,  $r=2M$", xy=(RH * math.cos(2.45), RH * math.sin(2.45)),
            xytext=(-2.42, 2.10), fontsize=9.4, color="#111111", zorder=45,
            arrowprops=dict(arrowstyle="-", color="#555555", lw=0.85, shrinkA=2, shrinkB=3))
ax.annotate("$r=M$", xy=(M * math.cos(0.72), M * math.sin(0.72)),
            xytext=(1.30, 1.34), fontsize=9.8, color="#111111", zorder=45,
            arrowprops=dict(arrowstyle="-", color="#555555", lw=0.85, shrinkA=2, shrinkB=3))
ax.plot([0], [0], marker="o", ms=4.4, color="#141414", zorder=46)
ax.annotate("singularity", xy=(0.02, -0.02), xytext=(0.30, -0.26), fontsize=9.4,
            color="#111111", zorder=45, va="center", ha="left",
            arrowprops=dict(arrowstyle="-", color="#555555", lw=0.85, shrinkA=2, shrinkB=3))

ax.annotate("", xy=(-M, -0.58), xytext=(M, -0.58), zorder=45,
            arrowprops=dict(arrowstyle="<->", color="#2b2b2b", lw=1.1))
ax.text(0, -0.74, "the two sheets can touch here", ha="center", va="top",
        fontsize=9.9, color="#141414", zorder=45)
ax.text(-2.44, -1.06, "one sixteenth of a turn to a band", ha="left", va="bottom",
        fontsize=8.3, color="#7a7a7a", zorder=45)

ax.set_xlim(-2.46, 2.30)
ax.set_ylim(-1.12, 2.22)
ax.set_aspect("equal")
ax.axis("off")

# ---------------------------------------------------------- right: the quantity
ax2 = fig.add_subplot(gs[1])
rr = np.linspace(1e-12, RH, 1200)
dd = 2 * np.pi - 4 * np.arcsin(np.sqrt(rr / RH))

ax2.fill_between(rr, np.pi, dd, where=(dd >= np.pi), color="#cbbce9", alpha=0.6, lw=0)
ax2.plot(rr, dd, color="#141414", lw=2.0, zorder=6)
ax2.axhline(np.pi, color="#8a2f2f", lw=1.4, ls=(0, (5, 3)), zorder=5)
ax2.axvline(M, color="#666666", lw=0.9, ls=(0, (2, 2.6)), zorder=4)

ax2.text(1.99, 3.26, "what the antipode asks for", ha="right", va="bottom",
         fontsize=8.9, color="#8a2f2f")
ax2.text(0.10, 3.40, "contact", ha="left", fontsize=9.7, color="#3b2a55")
ax2.text(1.45, 1.30, "no contact", ha="center", fontsize=9.7, color="#5c5c5c")

for k in range(NB):                                   # lock the axis to the same bands
    ax2.add_patch(plt.Rectangle((EDGES[k + 1], -0.40), EDGES[k] - EDGES[k + 1], 0.25,
                                facecolor=step_colour(k), edgecolor="white", lw=0.4,
                                clip_on=False, zorder=8))

ax2.set_xlim(0, RH)
ax2.set_ylim(-0.40, 2 * np.pi + 0.30)
ax2.set_yticks([0, np.pi, 2 * np.pi])
ax2.set_yticklabels(["$0$", r"$\pi$", r"$2\pi$"], fontsize=9.8)
ax2.set_xticks([0, M, RH])
ax2.set_xticklabels(["$0$", "$M$", "$2M$"], fontsize=9.8)
ax2.tick_params(axis="x", pad=11)
ax2.set_xlabel("radius", fontsize=9.8, labelpad=4)
ax2.set_ylabel("angle a causal curve can still sweep", fontsize=9.5)
ax2.tick_params(length=3)
for s_ in ("top", "right"):
    ax2.spines[s_].set_visible(False)
for s_ in ("left", "bottom"):
    ax2.spines[s_].set_color("#666666")

fig.text(0.025, 0.945, "The contact region, and the angle that decides it",
         fontsize=12.8, fontweight="bold", color="#111111")
fig.text(0.025, 0.882,
         "Colour is the angle itself. The spectrum starts exactly where contact does.",
         fontsize=9.7, color="#535353")

fig.savefig(os.path.join(OUT, "fig_companion_interior.pdf"))
fig.savefig(os.path.join(OUT, "fig_companion_interior.png"), dpi=200)
plt.close(fig)

# ------------------------------------------------------------------------ checks
P = print
P("fig_companion_interior: anchors")
P("   dphi_max(2M) = %.10f   want 0        %s" % (dphi_max(RH), abs(dphi_max(RH)) < 1e-12))
P("   dphi_max(M)  = %.10f   want pi       %s" % (dphi_max(M), abs(dphi_max(M) - math.pi) < 1e-12))
P("   dphi_max(0)  = %.10f   want 2pi      %s" % (dphi_max(0), abs(dphi_max(0) - 2 * math.pi) < 1e-12))

lo, hi = 1e-12, RH                      # solve dphi_max(r) = pi without using the closed form
for _ in range(200):
    mid = 0.5 * (lo + hi)
    if dphi_max(mid) > math.pi:
        lo = mid
    else:
        hi = mid
P("   crossing by bisection r = %.12f, want M = %.12f   %s"
  % (0.5 * (lo + hi), M, abs(0.5 * (lo + hi) - M) < 1e-10))

# A plant: the same figure drawn with 4 -> 3 must not cross pi at r = M.
wrong = 2 * math.pi - 3 * math.asin(math.sqrt(M / RH))
P("   plant: the 3-coefficient form misses pi at r=M: %s" % (abs(wrong - math.pi) > 0.1))
P("   wrote fig_companion_interior.pdf and .png")

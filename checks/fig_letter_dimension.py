"""How far round a causal curve can get, against how far the antipodal map asks it to go.

In D-dimensional Schwarzschild-Tangherlini the angle a causal curve can sweep on the two
interior crossings joining a point to its fold image is

    dphi_max(x) = (4/n) [ pi/2 - arcsin x^{n/2} ],    x = r/r_h,   n = D - 3,

rising to 2 pi/n at the singularity. The antipode asks for pi, flat in x and in D, because
that is fixed by the map being a free involution of a sphere. Contact happens wherever the
curve lies above the line, and only D = 4 ever does.

This replaces fig_letter_dimension.R. Two reasons. Its labels sat on the curves and on the
r = r_h/2 marker, the same fault Ben caught in the Penrose diagram and which no R script
here can check; and it was named fig_dimension_budget.R and its own comments priced the
antipodal map, which is a banned register. The layout is now measured by figaudit.Panel
and the script refuses to write a figure that collides.
"""
import math
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figaudit import Panel

OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                   "pub", "paper2")

ASKED = "#a8400f"
COLS = {4: "#1f4e79", 5: "#2e7d5b", 6: "#8a6d1f", 7: "#7a4a7a"}


def dphi(x, n):
    """Angle available on the two interior crossings, from depth x = r/r_h."""
    return (4.0 / n) * (math.pi / 2.0 - np.arcsin(np.asarray(x) ** (n / 2.0)))


fig, ax = plt.subplots(figsize=(7.0, 3.6))
fig.patch.set_facecolor("white")
p = Panel(fig, ax)

x = np.linspace(0, 1, 900)
ax.add_patch(plt.Rectangle((0, math.pi), 1, math.pi, facecolor="#f7ece5",
                           edgecolor="none", zorder=0))
p.line([0, 1], [math.pi, math.pi], "the line", color=ASKED, lw=2.0, zorder=3)
for D in (7, 6, 5, 4):
    p.line(x, dphi(x, D - 3), "D=%d curve" % D, color=COLS[D], lw=2.0, zorder=4,
           label="$D=%d$" % D)

ax.plot([0.5], [math.pi], marker="o", ms=6.5, color=COLS[4], zorder=6)
ax.plot([0.0], [math.pi], marker="o", ms=6.5, color=COLS[5], zorder=6)
p.line([0.5, 0.5], [0, math.pi], "r=rh/2 marker", color="#8a8a8a", lw=0.9,
       ls=(0, (2, 2.4)), zorder=2)

# A legend, not labels on the curves. Above the line only D = 4 exists and the upper right
# is empty, while below it four curves converge into the corner and every per-curve label
# landed on a neighbour: twelve collisions on the first attempt.
leg = ax.legend(loc="upper right", bbox_to_anchor=(0.995, 0.995), frameon=False,
                fontsize=9.2, handlelength=1.6, labelspacing=0.35, borderpad=0.2)
for txt, D in zip(leg.get_texts(), (7, 6, 5, 4)):
    txt.set_color(COLS[D])

p.text(0.985, 3.34, "what the antipodal map asks for", color=ASKED, fontsize=8.8,
       ha="right", va="bottom")
# The band between the line and the D = 4 curve is about 1.1 units tall at its widest and
# a two-line block does not fit in it. Above the D = 4 curve the whole shaded region is
# empty, so the text goes there instead.
p.text(0.44, 4.75, "$D=4$ clears it from $r=r_h/2$ inward:\nthe two sheets can touch",
       fontsize=8.6, color="#3b3b3b", ha="center", va="center", linespacing=1.4)
# Left of the r_h/2 marker and under the D = 7 curve is the only genuinely empty ground:
# right of it every curve sweeps down into the corner, and a box centred past about 0.33
# reaches the marker at 0.5. The first grid never got left of that and found nothing.
GRID = [(gx, gy) for gy in (0.34, 0.44, 0.54, 0.26, 0.64, 0.76)
        for gx in (0.30, 0.26, 0.32, 0.22, 0.34, 0.18, 0.74, 0.66)]
p.place(GRID, "every other dimension falls short\nat every depth",
        fontsize=8.6, color="#3b3b3b", ha="center", va="center", linespacing=1.4)
p.place([(gx, gy) for gy in (2.86, 2.70, 2.55, 2.40, 2.25) for gx in (0.16, 0.20, 0.24, 0.13)],
        "$D=5$ reaches it only here", fontsize=8.4, color=COLS[5],
        ha="center", va="center")

ax.set_xlim(0, 1)
ax.set_ylim(0, 2 * math.pi + 0.12)
ax.set_xticks([0, 0.25, 0.5, 0.75, 1.0])
ax.set_xticklabels(["0", "0.25", "0.5", "0.75", "1"], fontsize=9.4)
ax.set_yticks([0, math.pi / 2, math.pi, 3 * math.pi / 2, 2 * math.pi])
ax.set_yticklabels(["$0$", r"$\pi/2$", r"$\pi$", r"$3\pi/2$", r"$2\pi$"], fontsize=9.4)
ax.set_xlabel("depth  $r/r_h$   (0 = singularity, 1 = horizon)", fontsize=9.6, labelpad=5)
ax.set_ylabel("angle a causal curve can sweep", fontsize=9.6)
ax.tick_params(length=3)
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
for s in ("left", "bottom"):
    ax.spines[s].set_color("#666666")
fig.subplots_adjust(left=0.105, right=0.985, top=0.975, bottom=0.165)

P = print
P("fig_letter_dimension")
P("   plant: a label dropped on a curve is caught: %s" % p.selftest())
faults = p.audit()
if faults:
    P("   LAYOUT FAULTS (%d):" % len(faults))
    for f in faults:
        P("      " + f)
else:
    P("   layout clean: no text on a line, no text on text, nothing off canvas")

anchors = [("dphi(0.5, n=1) = pi", dphi(0.5, 1), math.pi),
           ("dphi(0, n=1) = 2pi", dphi(1e-300, 1), 2 * math.pi),
           ("dphi(0, n=2) = pi", dphi(1e-300, 2), math.pi),
           ("dphi(0, n=3) = 2pi/3", dphi(1e-300, 3), 2 * math.pi / 3)]
ok = True
for name, got, want in anchors:
    good = abs(float(got) - want) < 1e-9
    ok &= good
    P("   %-22s %.9f want %.9f  %s" % (name, got, want, good))
P("   plant: the D=5 curve never exceeds pi: %s"
  % bool(np.all(dphi(np.linspace(1e-12, 1, 500), 2) <= math.pi + 1e-12)))

if faults or not ok or not p.selftest():
    raise SystemExit("   figure NOT written: fix the faults above")
fig.savefig(os.path.join(OUT, "fig_letter_dimension.pdf"))
fig.savefig(os.path.join(OUT, "fig_letter_dimension.png"), dpi=220)
plt.close(fig)
P("   wrote fig_letter_dimension.pdf and .png")

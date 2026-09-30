"""Why a path from a point to its fold image crosses an interior exactly twice.

What does the work is the pair of shaded light cones. The future of any point in the
black-hole interior is a triangle capped by the singularity that reaches no horizon:
nothing leaves. The past of any point in the white-hole interior is capped the same way:
nothing re-enters. So a causal path from an image to its point meets the white-hole
interior once, an exterior, and the black-hole interior once, and no third crossing is
available to it. Visible, rather than claimed.

Coordinates. U and V are Kruskal, compactified by

    T = (arctan V + arctan U)/pi,      X = (arctan V - arctan U)/pi,

so the horizons UV = 0 are T = +-X, the singularities UV = 1 are the segments T = +-1/2
over |X| <= 1/2, and spatial infinity is X = +-1. Null rays run at forty-five degrees.
The fold's J:(U,V) -> (-U,-V) becomes (T,X) -> (-T,-X), reflection through the bifurcation
surface, so a point's image is diametrically opposite and the black-hole interior maps to
the white-hole interior.

The one easement, stated in the caption: exactly, r = M sits at T = 0.469 against a
singularity at 0.5, so a truthful drawing makes the contact region an unreadable sliver.
Constant-r curves are drawn at eased heights. The null structure is exact.

ON CHECKING THIS FIGURE. Two rounds of "the labels look fine to me" produced a figure with
text running over the diagram's own lines, and it took Ben reading the PDF to catch it.
So the layout is now measured, not judged: audit() renders the figure, takes every text
artist's bounding box and every drawn line's points in display coordinates, and reports
any text sitting on a line, any two texts sitting on each other, and any text outside the
canvas. A planted label proves the audit can fail. The script exits non-zero if anything
collides, so a bad layout cannot be written out and called finished.
"""
import math
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Polygon

OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                   "pub", "paper2")

INK, MUTED, EDGE = "#141414", "#6f6f6f", "#333333"
# One band per sheet, and the colours carry the physics rather than decorating it. The
# fold reverses time orientation, so the two contact regions are exact time-reverses of
# one another, and under time reversal redshift and blueshift exchange. Light escaping
# the black-hole side reaches an exterior redshifted; the same factor reversed makes
# light falling toward the white-hole side arrive blueshifted. A single mauve for both
# threw that away, which is what Ben noticed.
RED, BLUE = "#e9bfb6", "#bdcde6"          # black-hole side, white-hole side
# The points are ink, not red: red and blue now mean redshift and blueshift and nothing
# else, and a dark red dot on a pale red band was hard to see anyway.
PATH, PT, CONE = "#14406e", "#111111", "#8e8e8e"
BIGM, TMIN = 1.0, 0.30

LINES, TEXTS = [], []           # everything the audit measures


def _densify(xs, ys, per_unit=900):
    """Sample the whole line, not its corners. The hexagon is seven points, and a caption
    resting on one of its long edges was reported clean because nothing was sampled there."""
    xs, ys = np.asarray(xs, float), np.asarray(ys, float)
    if xs.size < 2:
        return xs, ys
    seg = np.hypot(np.diff(xs), np.diff(ys))
    total = float(seg.sum())
    if total <= 0:
        return xs, ys
    t = np.concatenate([[0.0], np.cumsum(seg)])
    u = np.linspace(0.0, total, max(xs.size, int(per_unit * total)))
    return np.interp(u, t, xs), np.interp(u, t, ys)


def line(ax, xs, ys, name, **kw):
    LINES.append((name,) + _densify(xs, ys))
    return ax.plot(xs, ys, **kw)


def label(ax, *a, **kw):
    t = ax.text(*a, **kw)
    TEXTS.append(t)
    return t


def note(ax, *a, **kw):
    t = ax.annotate(*a, **kw)
    TEXTS.append(t)
    return t


# ---------------------------------------------------------------------- the audit
def _text_rect(t, r, pad=(1.04, 1.14)):
    """Centre, width, height and angle of a text's real footprint, in display coords."""
    th = math.radians(t.get_rotation())
    keep = t.get_rotation()
    t.set_rotation(0)
    bb = t.get_window_extent(renderer=r)
    t.set_rotation(keep)
    # Rotation happens about the anchor, and the anchor is the box centre only when both
    # alignments are "center". Taking it as the centre regardless put the rectangle in the
    # wrong place for every right- or left-aligned label.
    ax_, ay = t.get_transform().transform(t.get_position())
    dx, dy = bb.x0 + bb.width / 2.0 - ax_, bb.y0 + bb.height / 2.0 - ay
    c, sn = math.cos(th), math.sin(th)
    cen = (ax_ + dx * c - dy * sn, ay + dx * sn + dy * c)
    return cen, bb.width * pad[0], bb.height * pad[1], th


def audit(planted=None):
    """Measure the layout. Returns a list of collisions; empty means clean."""
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    items = list(TEXTS) + ([planted] if planted is not None else [])
    boxes = [(t, t.get_window_extent(renderer=r).expanded(1.03, 1.10)) for t in items]
    bad = []
    for t, bb in boxes:
        txt = " ".join(t.get_text().split())[:30]
        cen, w, h, th = _text_rect(t, r)
        for name, xs, ys in LINES:
            pts = ax.transData.transform(np.column_stack([xs, ys]))
            # The true rotated rectangle, not the axis-aligned box matplotlib hands back.
            # For text at forty-five degrees the AABB is about (w+h)/sqrt(2) on a side,
            # far larger than the glyphs, and testing against it rejected every label that
            # ran along a diagonal. That is why the region names were not on the diagonals.
            du = pts - np.asarray(cen)
            u = du[:, 0] * math.cos(th) + du[:, 1] * math.sin(th)
            v = -du[:, 0] * math.sin(th) + du[:, 1] * math.cos(th)
            if np.any((np.abs(u) <= w / 2.0) & (np.abs(v) <= h / 2.0)):
                bad.append("text %-30s sits on %s" % ('"%s"' % txt, name))
        fb = fig.bbox
        if not (bb.x0 >= fb.x0 and bb.x1 <= fb.x1 and bb.y0 >= fb.y0 and bb.y1 <= fb.y1):
            bad.append("text %-30s runs off the canvas" % ('"%s"' % txt))
    for i in range(len(boxes)):
        for j in range(i + 1, len(boxes)):
            if boxes[i][1].overlaps(boxes[j][1]):
                bad.append('text "%s" overlaps "%s"'
                           % (" ".join(boxes[i][0].get_text().split())[:26],
                              " ".join(boxes[j][0].get_text().split())[:26]))
    return sorted(set(bad))


def place(ax, candidates, *a, **kw):
    """Put a text at the first candidate that collides with nothing.

    A candidate is (x, y) or (x, y, rotation), so a label can be offered one diagonal and
    then the other without the caller deciding in advance which has room.
    """
    t = ax.text(candidates[0][0], candidates[0][1], *a, **kw)
    TEXTS.append(t)
    key = '"%s"' % " ".join(t.get_text().split())[:30]
    for cand in candidates:
        t.set_position((cand[0], cand[1]))
        if len(cand) > 2:
            t.set_rotation(cand[2])
        if not [h for h in audit() if key in h]:
            return t
    return t


def along_diagonal(sign):
    """Candidates hugging an interior's own diagonal edges, nearest the edge first.

    Both edges are offered: the left one running up-left (rotation -45 above, +45 below)
    and the right one running up-right. Perpendicular clearance only has to beat half the
    label's height now that the audit measures the rotated rectangle instead of the
    axis-aligned box it sits in.
    """
    out = []
    for d in (0.045, 0.055, 0.065, 0.075, 0.090, 0.105):
        for m in (0.30, 0.26, 0.34, 0.22, 0.38, 0.18):
            out.append((-m, sign * (m + d), -45 if sign > 0 else 45))   # left edge
            out.append((m, sign * (m + d), 45 if sign > 0 else -45))    # right edge
    return out


def exact_rM_height(r=BIGM, m=BIGM):
    c = (1.0 - r / (2.0 * m)) * math.exp(r / (2.0 * m))
    return 2.0 * math.atan(math.sqrt(c)) / math.pi


def arc(sign, tmin=TMIN, n=260):
    x = np.linspace(-0.5, 0.5, n)
    return x, sign * (tmin + (0.5 - tmin) * (2.0 * x) ** 2)


def interior(x, t):
    return abs(t) > abs(x)


fig, ax = plt.subplots(figsize=(7.2, 4.0))
fig.patch.set_facecolor("white")

# ------------------------------------------------------------ contact bands, r < M
for sgn, col in ((+1, RED), (-1, BLUE)):
    x, t = arc(sgn)
    ax.fill_between(x, t, sgn * 0.5, color=col, alpha=0.85, lw=0, zorder=1)

# --------------------------------------- light cones: the reason for "and no more"
CONES = [(-0.18, 0.26, +1), (0.18, -0.26, -1)]
for cx0, ct0, sgn in CONES:
    reach = 0.5 - abs(ct0)
    ax.add_patch(Polygon([(cx0, ct0), (cx0 - reach, ct0 + sgn * reach),
                          (cx0 + reach, ct0 + sgn * reach)], closed=True,
                         facecolor=CONE, alpha=0.32, edgecolor="none", zorder=2))
    for s in (-1, +1):
        line(ax, [cx0, cx0 + s * reach], [ct0, ct0 + sgn * reach], "cone edge",
             color=CONE, lw=1.0, zorder=3)
    ax.plot([cx0], [ct0], marker="o", ms=3.0, color="#555555", zorder=4)

# ------------------------------------------------------------------- the diagram
B = [(1, 0), (0.5, 0.5), (-0.5, 0.5), (-1, 0), (-0.5, -0.5), (0.5, -0.5), (1, 0)]
line(ax, [p[0] for p in B], [p[1] for p in B], "boundary", color=EDGE, lw=1.5, zorder=6)
for a in (1, -1):
    line(ax, [-0.5, 0.5], [a * 0.5, a * -0.5], "horizon", color=EDGE, lw=1.1, zorder=5)
for sgn in (+1, -1):
    xz = np.linspace(-0.5, 0.5, 240)
    line(ax, xz, sgn * (0.5 + 0.011 * np.sin(48 * xz)), "singularity", color=INK,
         lw=2.2, zorder=7)
    x, t = arc(sgn)
    line(ax, x, t, "r=M", color=INK, lw=1.3, ls=(0, (5, 3)), zorder=7)

# ---------------------------------------------------------------- connecting path
XP, TP = 0.26, 0.46
knots = np.array([[-XP, -TP], [-0.04, -0.22], [0.12, 0.0], [0.20, 0.24], [XP, TP]])
u = np.linspace(0, 1, 420)
cx = np.interp(u, np.linspace(0, 1, len(knots)), knots[:, 0])
cy = np.interp(u, np.linspace(0, 1, len(knots)), knots[:, 1])
w = np.hanning(23); w /= w.sum()
cx = np.convolve(np.pad(cx, 11, mode="edge"), w, mode="same")[11:-11]
cy = np.convolve(np.pad(cy, 11, mode="edge"), w, mode="same")[11:-11]
line(ax, cx, cy, "path", color=PATH, lw=2.2, zorder=9, solid_capstyle="round")
for sx, sy in ((-XP, -TP), (XP, TP)):
    ax.plot([sx], [sy], marker="o", ms=6.4, color=PT, zorder=10,
            markeredgecolor="white", markeredgewidth=1.0)
label(ax, XP + 0.05, TP, "$x$", fontsize=12, color=PT, zorder=11, va="center")
label(ax, -XP - 0.05, -TP, r"$\Theta x$", fontsize=12, color=PT, zorder=11,
      va="center", ha="right")

inside = np.array([interior(a, b) for a, b in zip(cx, cy)])
runs, i = [], 0
while i < len(inside):
    if inside[i]:
        j = i
        while j < len(inside) and inside[j]:
            j += 1
        runs.append((i, j)); i = j
    else:
        i += 1
for k, (i0, i1) in enumerate(runs, 1):
    m = (i0 + i1) // 2
    ax.plot([cx[m]], [cy[m]], marker="o", ms=12.0, color=PATH, zorder=11)
    ax.text(cx[m], cy[m], str(k), color="white", fontsize=7.8, ha="center", va="center",
            zorder=12, fontweight="bold")          # sits on the path on purpose

# ------------------------------------------------------------------- annotation
label(ax, 0, 0.545, "singularity", ha="center", va="bottom", fontsize=9.2, color=INK)
label(ax, 0, -0.555, "singularity", ha="center", va="top", fontsize=9.2, color=INK)
# "the sheets touch" is gone from the figure. The band it named is bounded by the labelled
# r = M curve and the caption says what the shading is, and every free slot inside the band
# was taken by the cone, its label, the path or x. A label with nowhere to sit is a label
# that lands on something.
# r = M goes directly under the arc's lowest point, where the interior is empty. In the
# top-right corner the boundary, a horizon, the singularity and the arc all converge and a
# label there sits on all four of them.
# r = M moves off the axis and under the right of the arc, which frees the middle of the
# interior for the region name. With it in the centre there was no clean slot anywhere.
# Outside the diagram, beside the end of the arc. Inside, the arc and a horizon converge
# and the gap between them is thinner than the label.
place(ax, along_diagonal(+1), "black hole", fontsize=7.2, color=MUTED,
      ha="center", va="center", zorder=8)
place(ax, along_diagonal(-1), "white hole", fontsize=7.2, color=MUTED,
      ha="center", va="center", zorder=8)
place(ax, [(gx, gt) for gt in (0.245, 0.230, 0.260, 0.215) for gx in (0.0, 0.03, -0.03, 0.06)]
      + [(-0.60, 0.50), (0.62, 0.50)],
      "$r=M$", fontsize=9.0, color=INK, ha="center", va="center", zorder=11)
place(ax, [(0.34 + 0.02 * k, 0.24 - 0.012 * k) for k in range(9)],
      "$r=2M$", fontsize=8.6, color=MUTED, rotation=45, ha="center", va="center", zorder=8)
# Ben, 2026-09-30: naming the regions on their own diagonals saves the reader linking
# "future cone" back to a black hole named only in a caption across the diagram.
# The cone captions live in the widest part of each exterior, near T = 0; higher up the
# wedge narrows and a three-line block runs out through the boundary. They also name the
# two interiors, which saves two more labels that had nowhere to sit: "black hole" fitted
# only in the narrow tip of the triangle and "white hole" landed on the path.
# No leaders. A straight one reads as a null ray and a curved dotted one reads as broken,
# so each cone carries its own name and the caption beside it carries the consequence.
place(ax, [(-0.18 + dx, 0.26 + dt) for dt in (0.160, 0.150, 0.170, 0.145, 0.175)
           for dx in (0.0, 0.01, -0.01, 0.02)],
      "future cone", ha="center", va="center", fontsize=7.0, color="#3d3d3d", zorder=11, linespacing=1.35)
label(ax, -0.500, 0.000, "every future cone inside\nthe black hole ends on the\nsingularity, so nothing leaves", fontsize=8.3, color="#4a4a4a", zorder=11,
      va="center", ha="center", linespacing=1.4)
place(ax, [(0.18 + dx, -0.26 - dt) for dt in (0.160, 0.150, 0.170, 0.145, 0.175)
           for dx in (0.0, -0.01, 0.01, -0.02)],
      "past cone", ha="center", va="center", fontsize=7.0, color="#3d3d3d", zorder=11, linespacing=1.35)
label(ax, 0.500, 0.000, "every past cone inside the\nwhite hole starts on the\nsingularity, so nothing re-enters", fontsize=8.3, color="#4a4a4a", zorder=11,
      va="center", ha="center", linespacing=1.4)
# The "exterior" labels are gone. A wedge cannot hold a three-line caption and a label
# as well, and the captions already name both interiors while the figure caption says
# the path crosses an exterior. Four collisions came from insisting on them.

ax.set_xlim(-1.10, 1.10)
ax.set_ylim(-0.64, 0.64)
ax.set_aspect("equal")
ax.axis("off")
fig.subplots_adjust(left=0.005, right=0.995, top=0.995, bottom=0.005)




P = print
P("fig_letter_penrose")

probe = ax.text(0.0, 0.0, "PLANTED LABEL ON THE BIFURCATION POINT", ha="center",
                fontsize=9, color="none")
planted_hits = audit(planted=probe)
caught = any("PLANTED" in h for h in planted_hits)
probe.remove()
P("   plant: a label dropped on the horizons is caught: %s" % caught)

hits = audit()
if hits:
    P("   LAYOUT FAULTS (%d):" % len(hits))
    for h in hits:
        P("      " + h)
else:
    P("   layout clean: no text on a line, no text on text, nothing off canvas")

P("   r=M exactly at T=%.4f, singularity at 0.5000; drawn at %.2f, caption says eased"
  % (exact_rM_height(), TMIN))
slope = float(np.max(np.abs(np.diff(cx)) / np.maximum(np.abs(np.diff(cy)), 1e-15)))
P("   steepest |dX/dT| on the path: %.3f   causal everywhere: %s" % (slope, slope <= 1.0))
vis = sorted({("white " if b < 0 else "black ") + "interior" if interior(a, b)
              else "exterior" for a, b in zip(cx, cy)})
P("   regions visited: %s" % ", ".join(vis))
P("   exactly two interior crossings: %s (found %d)" % (len(runs) == 2, len(runs)))
allok = True
for cx0, ct0, sgn in CONES:
    reach = 0.5 - abs(ct0)
    ss = np.linspace(0, reach, 400)[1:]
    meets = any(abs(abs(cx0 + d * s) - abs(ct0 + sgn * s)) < 1e-4
                for d in (-1, 1) for s in ss)
    allok &= (abs(cx0 - reach) <= 0.5 + 1e-9 and abs(cx0 + reach) <= 0.5 + 1e-9
              and not meets)
P("   both cones end on a singularity and meet no horizon: %s" % allok)

if hits or not caught or not allok or slope > 1.0 or len(runs) != 2:
    raise SystemExit("   figure NOT written: fix the faults above")

fig.savefig(os.path.join(OUT, "fig_letter_penrose.pdf"))
fig.savefig(os.path.join(OUT, "fig_letter_penrose.png"), dpi=220)
plt.close(fig)
P("   wrote fig_letter_penrose.pdf and .png")

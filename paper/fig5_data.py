"""Paper 2 -- Figure 5: Sigma m_nu and the dark-energy stance against DESI DR2.

Hand-generated SVG, matching the house style of compute_figs.py (matplotlib is not
present in the project venv: numpy 2.5.3, scipy 1.18.1 checked, no matplotlib).
PNG is made from the SVG by headless google-chrome, same flags as make_pngs.sh.

Every number below is transcribed from a specific line of a working note in the author's
private development record (not included in this release); the note and the arXiv source
sit next to the number.  Nothing here is fit, guessed, or
carried over from pretrained knowledge of the literature.

LEFT PANEL -- Sigma m_nu (meV), 0-200:
  beyond/push/P2_NUR_DM/P2_REPORT.md            section 4, lines 344-383 (the two
                                                 masses, and Sec 4.2's DESI DR2 bounds,
                                                 "all bounds verified from the cited
                                                 abstracts")
  calc/tangents/numass_mirror/numass_mirror_output.txt   section 8 (lines 267-283),
                                                 the machine-readable Delta-chi^2(Sigma
                                                 m_nu) profile for LCDM, w0waCDM and the
                                                 two-sheet toy -- geometry only (BAO +
                                                 Pantheon+ + omega_cb + theta_* priors,
                                                 explicitly NO CMB lensing or primary-
                                                 amplitude information: same file, lines
                                                 99-105).

RIGHT PANEL -- the (w0, wa) CPL plane:
  working note TAU_NUMASS_20260910.md     lines 729-792, three DESI-DR2-BAO+CMB
                                                 best fits with quoted supernova samples,
                                                 fetched verbatim from their abstracts
                                                 from the cited abstracts.

No uncertainty is invented anywhere: every error bar plotted is a 1-sigma figure quoted
verbatim in one of the two source notes above; where a bound in the literature carries no
quoted sigma (the base DESI DR2 Sigma m_nu bounds), it is drawn as a point value with no
error bar, by design.
"""
import subprocess
import numpy as np

from pathlib import Path
OUT = str(Path(__file__).resolve().parent) + "/"

# =====================================================================================
# STYLE -- copied from compute_figs.py's SVG section so Figure 5 matches Figures 1-4.
# =====================================================================================
INK, MID, FAINT = "#15181d", "#5b6672", "#c8cdd4"
BLUE, RED, GREEN, ORANGE, PURPLE = "#1f5fa8", "#b03a2e", "#1e6b45", "#c8721a", "#6b3fa0"
SUB = lambda t: f"<tspan baseline-shift='sub' font-size='9'>{t}</tspan>"
SUP = lambda t: f"<tspan baseline-shift='super' font-size='9'>{t}</tspan>"

def svg_open(W, H, title=None, subtitle=()):
    s = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" '
         f'viewBox="0 0 {W} {H}" font-family="Helvetica,Arial,sans-serif">',
         f'<rect width="{W}" height="{H}" fill="#ffffff"/>']
    if title:
        s.append(f'<text x="26" y="30" font-size="15.5" fill="{INK}" font-weight="600">{title}</text>')
    for i, line in enumerate(subtitle):
        s.append(f'<text x="26" y="{50+i*16}" font-size="11.5" fill="{MID}">{line}</text>')
    return s

def txt(s, x, y, t, size=12, anchor="start", fill=INK, weight="normal", rot=None, style=""):
    tr = f' transform="rotate({rot} {x:.1f} {y:.1f})"' if rot is not None else ''
    st = f' font-style="{style}"' if style else ''
    s.append(f'<text x="{x:.1f}" y="{y:.1f}" font-size="{size}" text-anchor="{anchor}" '
             f'fill="{fill}" font-weight="{weight}"{st}{tr}>{t}</text>')

def footer(s, W, H, lines):
    for i, ln in enumerate(lines):
        txt(s, 26, H - 16*(len(lines)-i) + 4, ln, 10.3, fill=MID)

def vline(s, x, y0, y1, col, width=2.0, dash=None):
    da = f' stroke-dasharray="{dash}"' if dash else ""
    s.append(f'<line x1="{x:.1f}" y1="{y0:.1f}" x2="{x:.1f}" y2="{y1:.1f}" '
             f'stroke="{col}" stroke-width="{width}"{da}/>')

def hline(s, x0, x1, y, col, width=2.0, dash=None):
    da = f' stroke-dasharray="{dash}"' if dash else ""
    s.append(f'<line x1="{x0:.1f}" y1="{y:.1f}" x2="{x1:.1f}" y2="{y:.1f}" '
             f'stroke="{col}" stroke-width="{width}"{da}/>')

def polyline(s, pts, col, width=2.2, dash=None):
    da = f' stroke-dasharray="{dash}"' if dash else ""
    p = " ".join(f"{x:.2f},{y:.2f}" for x, y in pts)
    s.append(f'<polyline points="{p}" fill="none" stroke="{col}" stroke-width="{width}"{da}/>')

SMNU = "&#931;m" + SUB("&#957;")                # Sigma m_nu
CHI2 = "&#916;&#967;" + SUP("2")                # Delta chi^2
W0   = "w" + SUB("0")
WA   = "w" + SUB("a")

# =====================================================================================
# LEFT-PANEL DATA
# =====================================================================================
# The two framework values (BFT seesaw, one vanishing Yukawa row => one massless
# neutrino => Sigma m_nu at the oscillation minimum for that ordering).
# NuFIT 6.0 Table 1, IC24+SK central values --
#   "| NO, m_1 = 0 | (0, 8.65, 50.13) | 58.8 meV  | 1.50-3.70 meV |"
#   "| IO, m_3 = 0 | (50.10,50.85,0)  | 100.9 meV | 18.2-48.2 meV |"
NH_SUM = 58.8     # meV -- normal ordering, m1 = 0 -- the framework's live prediction
IH_SUM = 100.9     # meV -- inverted ordering, m3 = 0 -- excluded (see bound below)

# DESI DR2 bounds, all Elbers et al. 2503.14744, quoted in P2_REPORT.md Sec 4.2:
#   P2_REPORT.md:362-364  "DESI DR2 BAO + CMB, LCDM, Sigma m_nu >= 0 prior:
#                          Sigma m_nu < 0.0642 eV (95%)"                    -> 64.2 meV
#   P2_REPORT.md:365-366  "Feldman-Cousins boundary-corrected:
#                          Sigma m_nu < 0.053 eV (95%)"                     -> 53   meV
#   P2_REPORT.md:367      "w0waCDM: Sigma m_nu < 0.163 eV (95%)"            -> 163  meV
# (eV -> meV is a factor of 1000; the source states the eV figure, reproduced here as
# given, only rescaled to match the figure's meV axis.)
LCDM_UL  = 64.2   # meV, 95%, DESI DR2 BAO+CMB, LCDM     [2503.14744]
FC_UL    = 53.0   # meV, 95%, Feldman-Cousins correction  [2503.14744]
W0WA_UL  = 163.0  # meV, 95%, DESI DR2 BAO+CMB, w0waCDM  [2503.14744]

# Delta-chi^2(Sigma m_nu), geometry only, relative to EACH model's own minimum (each of
# the three models' own chi2_min occurs at Sigma m_nu = 0, so the tabulated "D_*"
# columns already are Delta-chi^2 relative to that model's minimum -- verified against
# numass_mirror_output.txt section 3, "chi2_min = ... at Sum m_nu = 0.0 meV" for all
# three).  Table transcribed verbatim from numass_mirror_output.txt:267-283:
#
#   Sum_meV  D_LCDM   D_w0wa   D_toy+   D_toy+-
SUM_GRID = np.array([0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 58.8, 70.0, 85.0,
                      100.0, 120.0, 150.0, 180.0, 220.0, 260.0, 300.0])
D_LCDM = np.array([0.000, 0.212, 0.642, 1.168, 1.760, 2.405, 3.012, 3.834, 5.014,
                    6.281, 8.098, 11.086, 14.373, 19.195, 24.490, 30.228])
D_W0WA = np.array([0.000, 0.096, 0.295, 0.541, 0.821, 1.129, 1.420, 1.817, 2.393,
                    3.018, 3.923, 5.436, 7.131, 9.667, 12.512, 15.656])
# "D_toy+" (Omega_f >= 0) and "D_toy+-" (Omega_f any sign) are identical at every grid
# point in the source table -- transcribed once, as D_TOY, per numass_mirror_output.txt:
# "toy2side : ... Delta chi^2 at Sum m_nu = 58.8 meV vs this model's own best = +3.645"
# matching the Of>=0 row exactly.
D_TOY = np.array([0.000, 0.294, 0.835, 1.476, 2.182, 2.941, 3.645, 4.588, 5.925,
                   7.342, 9.351, 12.602, 16.125, 21.221, 26.743, 32.664])
DCHI2_95 = 2.71   # numass_mirror_output.txt:80, "95% profile-likelihood bound (Dchi2 = 2.71 ...)"

# =====================================================================================
# RIGHT-PANEL DATA -- the (w0, wa) CPL plane
# =====================================================================================
# The framework's own dark-energy content is an exact cosmological constant: w0=-1,
# wa=0, no phantom crossing (PAPER2_v3.md sec 5.6). This is a definition of the
# framework's stance, not a data point.
LCDM_W0, LCDM_WA = -1.0, 0.0

# Three published DESI-DR2-BAO + CMB best fits, each combined with a different
# supernova sample, fetched verbatim (marked [C]) in
# working note TAU_NUMASS_20260910.md.  1-sigma uncertainties are exactly as
# quoted; asymmetric ones are kept asymmetric rather than symmetrised.
DESI_POINTS = [
    # label,                 arXiv id,     w0,     w0_lo,  w0_hi,  wa,    wa_lo, wa_hi
    # TAU_NUMASS_20260910.md:734-738 (arXiv:2511.07517, Popovic et al., DES-Dovekie
    # recalibration + CMB[Planck,ACT,SPT] + DESI DR2, flat w0waCDM):
    # "we find w0 = -0.803 +/- 0.054, wa = -0.72 +/- 0.21."
    ("DESI DR2+CMB+DES-Dovekie", "2511.07517", -0.803, 0.054, 0.054, -0.72, 0.21, 0.21),
    # TAU_NUMASS_20260910.md:758-760 (arXiv:2601.19424, Hoyt et al., host-mass-corrected
    # Union3.1 + "BAO ... and CMB exactly as done by DESI DR2"):
    # "we find w0=-0.719+/-0.084, wa=-0.95^{+0.29}_{-0.26}."
    ("DESI DR2+CMB+Union3.1(corr.)", "2601.19424", -0.719, 0.084, 0.084, -0.95, 0.26, 0.29),
    # TAU_NUMASS_20260910.md:786-789 (arXiv:2609.05053, "Supernovae Unite",
    # Pantheon+ & DES-SN5YR combined, + BAO + CMB):
    # "(Omega_m,w0,wa) = (0.305+/-0.004, -0.861^{+0.044}_{-0.042}, -0.60^{+0.17}_{-0.19})"
    ("DESI DR2+CMB+Unite(Pantheon+/DES-SN5YR)", "2609.05053", -0.861, 0.042, 0.044, -0.60, 0.19, 0.17),
]

# =====================================================================================
# FIGURE 5
# =====================================================================================
W, H = 1100, 610
s = svg_open(W, H,
    "Figure 5 &#8212; " + SMNU + " and the dark-energy stance against DESI DR2",
    ["Left: the seesaw prediction and the oscillation floor against DESI DR2 " + SMNU +
     " bounds [2503.14744], with this work's own geometry-only " + CHI2 + "(" + SMNU +
     ") reproduction (no CMB lensing).",
     "Right: the framework's exact-&#923; stance (no phantom crossing) against three "
     "published DESI DR2 + CMB (w&#8320;,w&#8336;) fits with their quoted 1&#963; errors."])

MT, PH = 132, 340

# ------------------------------------------------------------------ LEFT PANEL
ML, PW = 84, 436
XMIN, XMAX = 0.0, 200.0          # meV
YMIN, YMAX = 0.0, 10.0           # Delta chi^2, relative to each model's own minimum
XL = lambda v: ML + (v - XMIN) / (XMAX - XMIN) * PW
YL = lambda v: MT + PH - (v - YMIN) / (YMAX - YMIN) * PH

# oscillation-floor band: Sigma m_nu < 58.8 meV is kinematically impossible for the
# normal ordering (m1 >= 0 and the measured splittings fix the rest) -- P2_REPORT.md's
# own reading of this, verbatim at :370-371: an independent 53 meV bound sits "below
# the minimum of 59 meV set by the normal ordering" [2507.12401].
s.append(f'<rect x="{XL(XMIN):.1f}" y="{MT}" width="{XL(NH_SUM)-XL(XMIN):.1f}" height="{PH}" '
         f'fill="{FAINT}" opacity="0.55"/>')
s.append(f'<rect x="{ML}" y="{MT}" width="{PW}" height="{PH}" fill="none" stroke="{INK}"/>')

for t in (0, 50, 100, 150, 200):
    s.append(f'<line x1="{XL(t):.1f}" y1="{MT+PH}" x2="{XL(t):.1f}" y2="{MT+PH+5}" stroke="{INK}"/>')
    txt(s, XL(t), MT+PH+18, str(t), 10.5, "middle")
txt(s, ML+PW/2, MT+PH+38, SMNU + "   (meV)", 12.5, "middle")

for t in (0, 2, 4, 6, 8, 10):
    s.append(f'<line x1="{ML-5}" y1="{YL(t):.1f}" x2="{ML}" y2="{YL(t):.1f}" stroke="{INK}"/>')
    txt(s, ML-9, YL(t)+4, str(t), 10.5, "end")
txt(s, 22, MT+PH/2, CHI2 + "  (each model, rel. to its own minimum)", 11.5, "middle", rot=-90)

# Delta chi^2 = 2.71 reference line (drawn now for the axis; redrawn on top of the
# legend boxes below so it stays visible where it passes behind them, and labelled
# there as legend-box row 6 rather than as a floating label -- every open spot on this
# axis collides with either a vertical line or a curve at this data density).
hline(s, ML, ML+PW, YL(DCHI2_95), MID, 1.6, "5 4")

# the two framework predictions
vline(s, XL(NH_SUM), MT, MT+PH, BLUE, 2.6)
vline(s, XL(IH_SUM), MT, MT+PH, RED, 2.6)
s.append(f'<line x1="{XL(IH_SUM)-6:.1f}" y1="{MT+7}" x2="{XL(IH_SUM)+6:.1f}" y2="{MT+19}" stroke="{RED}" stroke-width="2"/>')
s.append(f'<line x1="{XL(IH_SUM)+6:.1f}" y1="{MT+7}" x2="{XL(IH_SUM)-6:.1f}" y2="{MT+19}" stroke="{RED}" stroke-width="2"/>')

# the three DESI DR2 bounds (point values; no uncertainty is quoted for these in the
# source, so none is drawn)
vline(s, XL(LCDM_UL), MT, MT+PH, BLUE, 2.0, "6 4")
vline(s, XL(FC_UL), MT, MT+PH, ORANGE, 2.0, "2 3")
vline(s, XL(W0WA_UL), MT, MT+PH, PURPLE, 2.0, "9 3 2 3")

# Delta chi^2(Sigma m_nu) curves -- piecewise-linear interpolation of the tabulated grid
# above (numass_mirror_output.txt:267-283); clipped to the panel once a curve exceeds
# Delta chi^2 = 10 (all three are monotonically increasing over this range).
xs = np.linspace(XMIN, XMAX, 401)
for ys, col in ((np.interp(xs, SUM_GRID, D_LCDM), BLUE),
                (np.interp(xs, SUM_GRID, D_W0WA), PURPLE),
                (np.interp(xs, SUM_GRID, D_TOY), GREEN)):
    mask = ys <= YMAX
    pts = list(zip(XL(xs[mask]), YL(ys[mask])))
    polyline(s, pts, col, 2.6)

# legend box 1: the vertical reference lines, plus the Delta-chi^2=2.71 threshold
lx, ly, lw, lh = ML+8, MT+8, 258, 148
s.append(f'<rect x="{lx}" y="{ly}" width="{lw}" height="{lh}" fill="#ffffff" opacity="0.93" stroke="{FAINT}"/>')
rows1 = [("NH, m&#8321;=0 (exact-rule sector): 58.8 meV", BLUE, None),
         ("IH, m&#8323;=0 (exact-rule sector): 100.9 meV", RED, None),
         ("DESI DR2+CMB &#923;CDM 95% UL: 64.2 meV", BLUE, "6 4"),
         ("Feldman&#8211;Cousins 95%: 53 meV", ORANGE, "2 3"),
         ("DESI DR2+CMB w&#8320;w&#8336;CDM 95% UL: 163 meV", PURPLE, "9 3 2 3"),
         (CHI2 + " = 2.71 (95%, profile likelihood)", MID, "5 4")]
for i, (lab, col, dash) in enumerate(rows1):
    yy = ly + 16 + i*20
    da = f' stroke-dasharray="{dash}"' if dash else ""
    s.append(f'<line x1="{lx+8}" y1="{yy}" x2="{lx+30}" y2="{yy}" stroke="{col}" stroke-width="2.4"{da}/>')
    txt(s, lx+38, yy+4, lab, 9.6, fill=INK)
txt(s, lx+8, ly+lh-6, "shaded: below the m&#8321;=0 floor (kinematically forbidden, NH)", 8.8, fill=MID)

# legend box 2: the Delta chi^2 curves -- right-aligned so it always sits inside the
# panel border regardless of PW (was overflowing when pinned to a data x-coordinate).
lw2, lh2 = 230, 86
lx2, ly2 = ML + PW - lw2 - 8, MT + PH - 96
s.append(f'<rect x="{lx2:.1f}" y="{ly2}" width="{lw2}" height="{lh2}" fill="#ffffff" opacity="0.93" stroke="{FAINT}"/>')
rows2 = [(CHI2 + "(" + SMNU + "): &#923;CDM", BLUE),
         (CHI2 + "(" + SMNU + "): w&#8320;w&#8336;CDM", PURPLE),
         (CHI2 + "(" + SMNU + "): two-sheet toy", GREEN)]
for i, (lab, col) in enumerate(rows2):
    yy = ly2 + 16 + i*20
    s.append(f'<line x1="{lx2+8:.1f}" y1="{yy}" x2="{lx2+30:.1f}" y2="{yy}" stroke="{col}" stroke-width="2.6"/>')
    txt(s, lx2+38, yy+4, lab, 9.6, fill=INK)
txt(s, lx2+8, ly2+86-6, "geometry only, no CMB lensing (weaker than bounds at left)", 8.6, fill=MID)

# redraw the Delta-chi^2=2.71 line on top of both legend boxes so it reads as one
# continuous guide across the full panel width, not truncated where it passes under them
hline(s, ML, ML+PW, YL(DCHI2_95), MID, 1.6, "5 4")

txt(s, ML+10, MT+PH-12, "(a)", 12, fill=INK, weight="600")

# ------------------------------------------------------------------ RIGHT PANEL
ML2, PW2 = 636, 424
W0MIN, W0MAX = -1.08, -0.60
WAMIN, WAMAX = -1.30, 0.25
XR = lambda v: ML2 + (v - W0MIN) / (W0MAX - W0MIN) * PW2
YR = lambda v: MT + PH - (v - WAMIN) / (WAMAX - WAMIN) * PH

s.append(f'<rect x="{ML2}" y="{MT}" width="{PW2}" height="{PH}" fill="none" stroke="{INK}"/>')
for t in (-1.0, -0.9, -0.8, -0.7, -0.6):
    s.append(f'<line x1="{XR(t):.1f}" y1="{MT+PH}" x2="{XR(t):.1f}" y2="{MT+PH+5}" stroke="{INK}"/>')
    txt(s, XR(t), MT+PH+18, f"{t:.1f}".replace("-", "&#8722;"), 10.5, "middle")
txt(s, ML2+PW2/2, MT+PH+38, W0, 12.5, "middle")
for t in (0.2, 0.0, -0.2, -0.4, -0.6, -0.8, -1.0, -1.2):
    s.append(f'<line x1="{ML2-5}" y1="{YR(t):.1f}" x2="{ML2}" y2="{YR(t):.1f}" stroke="{INK}"/>')
    lab = f"{t:+.1f}".replace("+", "").replace("-", "&#8722;")
    txt(s, ML2-9, YR(t)+4, lab, 10.5, "end")
txt(s, ML2-46, MT+PH/2, WA, 12.5, "middle", rot=-90)

# phantom-divide reference line, w0 = -1 (today)
vline(s, XR(-1.0), MT, MT+PH, MID, 1.6, "5 4")
txt(s, XR(-1.0)+6, MT+14, "w = &#8722;1 (phantom divide, today)", 9.6, fill=MID)

# three published DESI DR2 + CMB (w0, wa) fits, with quoted 1-sigma error bars
for i, (label, arxiv, w0, w0lo, w0hi, wa, walo, wahi) in enumerate(DESI_POINTS):
    cx, cy = XR(w0), YR(wa)
    s.append(f'<line x1="{XR(w0-w0lo):.1f}" y1="{cy:.1f}" x2="{XR(w0+w0hi):.1f}" y2="{cy:.1f}" '
             f'stroke="{INK}" stroke-width="1.6"/>')
    s.append(f'<line x1="{cx:.1f}" y1="{YR(wa-walo):.1f}" x2="{cx:.1f}" y2="{YR(wa+wahi):.1f}" '
             f'stroke="{INK}" stroke-width="1.6"/>')
    s.append(f'<circle cx="{cx:.1f}" cy="{cy:.1f}" r="4.2" fill="{INK}"/>')
    dy = -10 if i != 1 else 16
    txt(s, cx+8, cy+dy, f"{i+1}", 10.5, fill=INK, weight="600")

# the LCDM point -- the framework's stance
s.append(f'<circle cx="{XR(LCDM_W0):.1f}" cy="{YR(LCDM_WA):.1f}" r="6.4" fill="#ffffff" stroke="{BLUE}" stroke-width="2.6"/>')
s.append(f'<circle cx="{XR(LCDM_W0):.1f}" cy="{YR(LCDM_WA):.1f}" r="2.2" fill="{BLUE}"/>')
txt(s, XR(LCDM_W0)-8, YR(LCDM_WA)-10, "&#923;CDM (the framework&#8217;s stance)", 10, "end", BLUE, "600")

# legend: the three numbered DESI points. Placed in the panel's top-right corner,
# which the point cloud, its error bars, the toy curve and the phantom-divide label
# all leave clear (checked against each element's plotted pixel extent); short tags
# here, full sample names in fig5_caption.txt and in DESI_POINTS above.
SHORT_TAGS = ["Dovekie", "Union3.1(corr.)", "Unite"]
lx3, ly3, lw3, lh3 = ML2+PW2-186, MT+8, 178, 78
s.append(f'<rect x="{lx3}" y="{ly3}" width="{lw3}" height="{lh3}" fill="#ffffff" opacity="0.93" stroke="{FAINT}"/>')
for i, (label, arxiv, *_r) in enumerate(DESI_POINTS):
    yy = ly3 + 15 + i*22
    s.append(f'<circle cx="{lx3+10}" cy="{yy-4}" r="4.2" fill="{INK}"/>')
    txt(s, lx3+22, yy, f"{i+1}", 9.6, fill=INK, weight="600")
    txt(s, lx3+36, yy, f"{SHORT_TAGS[i]}  [{arxiv}]", 8.8, fill=INK)

txt(s, ML2+10, MT+PH-12, "(b)", 12, fill=INK, weight="600")

# ------------------------------------------------------------------ FOOTER
footer(s, W, H, [
 "Sources: " + SMNU + " values and DESI DR2 bounds -- beyond/push/P2_NUR_DM/P2_REPORT.md "
 "&#167;4.2 (bounds all from Elbers et al., arXiv:2503.14744). " + CHI2 + "(" + SMNU + ") profiles --",
 "calc/tangents/numass_mirror/numass_mirror_output.txt &#167;8 (geometry-only: BAO + Pantheon+ + "
 "&#969;" + SUB("cb") + " + &#952;" + SUB("*") + " priors, no CMB lensing or primary-amplitude term).",
 "(w&#8320;,w&#8336;) points -- working note TAU_NUMASS_20260910.md, arXiv:2511.07517, "
 "2601.19424 and 2609.05053."])

s.append("</svg>")
open(OUT + "fig5_data.svg", "w").write("\n".join(s))
print("wrote fig5_data.svg")

# =====================================================================================
# RASTERISE -- same headless-chrome recipe as make_pngs.sh (matplotlib not available).
# =====================================================================================
cmd = ["google-chrome", "--headless", "--disable-gpu", "--no-sandbox", "--hide-scrollbars",
       "--force-device-scale-factor=2", f"--window-size={W},{H}",
       f"--screenshot={OUT}fig5_data.png", f"file://{OUT}fig5_data.svg"]
r = subprocess.run(cmd, capture_output=True, text=True)
if r.returncode != 0:
    print("google-chrome stderr:", r.stderr[-2000:])
else:
    import os
    print(f"wrote fig5_data.png  {os.path.getsize(OUT+'fig5_data.png')} bytes")

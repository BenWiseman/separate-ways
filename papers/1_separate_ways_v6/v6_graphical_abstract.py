#!/usr/bin/env python3
"""Graphical abstract for Separate Ways v6 (Nuclear Physics B).

The new result leads: a line of widely used graviton-exchange rates counts the plasma's fermions
twice, which makes them 1.41 times too large. A schematic on the left gives the setting (a universe
whose far side of the bang is its CPT image, with dark matter made in the hot plasma), the corrected
rates sit in the middle, and three cards on the right give what follows for the dark-matter mass.
Colours are the house palette of v6_style.py, chosen to stay distinct in the common forms of colour
blindness and in greyscale; every colour contrast also carries a label or a hatching.

It adds no number. Every value is parsed from the printed output of the script that makes it, the
same lines v6_numbers_check.py checks: the rate ratios from v6_grav_rate.py, the 8.0e13 GeV plasma
temperature from v6_bath.py, the 491.6 PeV ceiling and its 245.8 PeV neutrino line from
v6_abundance.py, and the neutrino-mass sum from v6_neutrino.py. The 480 PeV in the setting is Boyle,
Finn and Turok's own figure, as the abstract quotes it.

The page is HTML with inline SVG, 10 x 4 inches, rendered by headless Chrome to a PDF and to a
3000 x 1200 pixel PNG (Elsevier asks for at least 1328 x 531 pixels).
"""
import os, re, subprocess, sys

HERE = os.path.dirname(os.path.abspath(__file__))
CALC = os.path.join(HERE, "..", "..", "checks", "calc")
CHROME = os.environ.get("CHROME", "google-chrome")
STEM = os.path.join(HERE, "graphical_abstract_NPB_v6")

def run(script):
    r = subprocess.run([sys.executable, os.path.join(CALC, script)], capture_output=True, text=True, timeout=3600)
    if r.returncode != 0: raise SystemExit("%s exited %d" % (script, r.returncode))
    return r.stdout

def grab(pattern, text, what):
    m = re.search(pattern, text)
    if not m: raise SystemExit("no %s in the script's output" % what)
    return m

def numbers():
    g = run("v6_grav_rate.py")
    ratios = [(lab, float(grab(re.escape(key) + r"\s+r = .*?with 45 Dirac [0-9.e+-]+, ([0-9.]+) times too high", g, key).group(1)))
              for lab, key in (("scalar", "real scalar, xi = 0"), ("fermion", "Dirac fermion"), ("vector", "massless vector"))]
    T_star = float(grab(r"T\* \(alone matches the abundance at M_c\) = ([0-9.e+]+) GeV", run("v6_bath.py"), "T*").group(1))
    ab = run("v6_abundance.py")
    M_ceiling = float(grab(r"M_\* = ([0-9.]+) PeV", ab, "the ceiling").group(1))
    line = float(grab(r"line at M_\*/2 = ([0-9.]+) PeV", ab, "the line").group(1))
    nu = float(grab(r"normal ordering, m1 = 0:\s+Sum = ([0-9.]+)", run("v6_neutrino.py"), "the neutrino-mass sum").group(1))
    return ratios, T_star, M_ceiling, line, nu

# the house palette of v6_style.py: blue for corrected values and our universe, deep purple for the
# published count and the CPT image, orange for what to look for
sys.path.insert(0, HERE)
import v6_style as st
P_DEEP, P, P_MID, P_LIGHT, P_PALE = st.PURPLE_D, st.PURPLE, st.PURPLE_M, st.PURPLE_L, st.PURPLE_P
B_DEEP, B, B_MID, B_LIGHT, B_PALE = st.BLUE_D, st.BLUE, st.BLUE_M, st.BLUE_L, st.BLUE_P
O_DEEP, O, O_PALE = st.ORANGE_D, st.ORANGE, st.ORANGE_P
INK, GREY = "#1e1b24", "#5b5866"

def bars_svg(ratios, toc=False):
    """corrected (blue, = 1) against the published count (purple, hatched) for each kind of dark matter"""
    W, H, base, top = (330, 150, 118, 26) if toc else (340, 160, 124, 30)
    scale = (base - top)/1.6
    out = ['<svg width="%d" height="%d" viewBox="0 0 %d %d">' % (W, H, W, H),
           '<defs><pattern id="hatch" width="6" height="6" patternUnits="userSpaceOnUse" patternTransform="rotate(45)">'
           '<rect width="6" height="6" fill="white"/><line x1="0" y1="0" x2="0" y2="6" stroke="%s" stroke-width="2.4"/>'
           '</pattern></defs>' % P,
           '<line x1="18" y1="%d" x2="%d" y2="%d" stroke="%s" stroke-width="1"/>' % (base, W - 6, base, GREY),
           '<line x1="18" y1="%.1f" x2="%d" y2="%.1f" stroke="%s" stroke-width="0.8" stroke-dasharray="3 3"/>'
           % (base - scale, W - 6, base - scale, B_MID)]
    bw, gap = 32, 108
    for i, (lab, r) in enumerate(ratios):
        x0 = 40 + i*gap
        h1, h2 = scale, scale*r
        out.append('<rect x="%d" y="%.1f" width="%d" height="%.1f" fill="%s" rx="2"/>' % (x0, base - h1, bw, h1, B))
        out.append('<rect x="%d" y="%.1f" width="%d" height="%.1f" fill="white" stroke="%s" stroke-width="2.2" rx="2"/>'
                   % (x0 + bw + 5, base - h2 + 1, bw - 2, h2 - 1, P))
        out.append('<text x="%d" y="%.1f" text-anchor="middle" class="ratio">&#215;%.2f</text>' % (x0 + bw + 4 + bw/2, base - h2 - 6, r))
        out.append('<text x="%d" y="%d" text-anchor="middle" class="tick">%s</text>' % (x0 + bw + 2, base + 16, lab))
        if not toc:
            out.append('<text x="%d" y="%d" text-anchor="middle" class="tick">dark matter</text>' % (x0 + bw + 2, base + 28))
    out.append('</svg>')
    return "\n".join(out)

def hourglass_svg(up=None, down=None, side_px=12, width=236):
    """our side of the bang (blue) and its CPT image (purple), meeting at the bang, with the hot plasma there"""
    up = up or (B_PALE, B); down = down or (P_PALE, P)
    W, H, cx, cy = width, 232, 106, 116
    dots = [(-14, -6, B_MID), (10, -9, P_MID), (-5, 5, P_MID), (16, 4, B_MID), (-22, 2, P_LIGHT), (24, -2, B_LIGHT),
            (2, -14, B_LIGHT), (-4, 13, P_LIGHT), (7, 11, B_MID), (-12, -14, P_MID)]
    out = ['<svg width="%d" height="%d" viewBox="0 0 %d %d">' % (W, H, W, H), '<defs>',
           '<linearGradient id="gUp" x1="0" y1="1" x2="0" y2="0"><stop offset="0" stop-color="%s"/><stop offset="1" stop-color="%s"/></linearGradient>' % up,
           '<linearGradient id="gDn" x1="0" y1="0" x2="0" y2="1"><stop offset="0" stop-color="%s"/><stop offset="1" stop-color="%s"/></linearGradient>' % down,
           '<radialGradient id="glow"><stop offset="0" stop-color="#ffffff"/><stop offset="0.55" stop-color="#fff6d8"/><stop offset="1" stop-color="#fff6d8" stop-opacity="0"/></radialGradient>',
           '<marker id="ah" markerWidth="8" markerHeight="8" refX="6" refY="4" orient="auto"><path d="M0,0 L8,4 L0,8 z" fill="%s"/></marker>' % GREY,
           '</defs>',
           '<polygon points="%d,%d %d,%d %d,%d" fill="url(#gUp)"/>' % (cx, cy, cx - 88, 6, cx + 88, 6),
           '<polygon points="%d,%d %d,%d %d,%d" fill="url(#gDn)"/>' % (cx, cy, cx - 88, H - 6, cx + 88, H - 6),
           '<circle cx="%d" cy="%d" r="34" fill="url(#glow)"/>' % (cx, cy)]
    for dx, dy, c in dots:
        out.append('<circle cx="%d" cy="%d" r="3.2" fill="%s"/>' % (cx + dx, cy + dy, c))
    out += ['<text x="%d" y="%d" text-anchor="middle" class="side" fill="white" style="font-size:%dpx">our universe</text>' % (cx, 15 + side_px, side_px),
            '<text x="%d" y="%d" text-anchor="middle" class="side" fill="white" style="font-size:%dpx">its CPT image</text>' % (cx, H - 16, side_px),
            '<text x="%d" y="%d" text-anchor="end" class="bang">the bang</text>' % (cx - 40, cy + 4),
            '<path d="M %d,%d C %d,%d %d,%d %d,%d" fill="none" stroke="%s" stroke-width="1.4" marker-end="url(#ah)" marker-start="url(#ah)"/>'
            % (cx + 58, 50, cx + 114, 86, cx + 114, 146, cx + 58, H - 50, GREY),
            '<text x="%d" y="%d" class="cpt">CPT</text>' % (cx + 108, cy + 4)]
    out.append('</svg>')
    return "\n".join(out)

def page(ratios, T_star, M_ceiling, line, nu):
    rf = dict(ratios)
    mant, ex = ("%.1e" % T_star).split("e+")
    return """<!doctype html>
<html><head><meta charset="utf-8"><title>Separate Ways v6 graphical abstract</title>
<style>
@page { size: 10in 4in; margin: 0; }
html, body { margin: 0; padding: 0; }
body { width: 960px; height: 384px; position: relative; overflow: hidden; background: #ffffff;
       font-family: Arial, "Liberation Sans", Helvetica, sans-serif; color: %(INK)s; }
.abs { position: absolute; }
.head { left: 28px; top: 12px; width: 904px; font-size: 25px; font-weight: 700; letter-spacing: -0.2px; }
.head .em { color: %(P)s; }
.sub { left: 28px; top: 47px; width: 904px; font-size: 13px; color: %(GREY)s; }
.big { left: 296px; top: 80px; font-size: 64px; font-weight: 700; color: %(P)s; line-height: 1; letter-spacing: -1px; }
.bigcap { left: 474px; top: 82px; width: 172px; font-size: 12px; line-height: 1.28; }
.bigcap b { color: %(B_DEEP)s; }
.tag { font-size: 10.5px; color: %(GREY)s; }
.side { font-size: 12px; font-weight: 700; }
.bang { font-size: 11.5px; font-weight: 700; fill: %(INK)s; }
.cpt { font-size: 11px; font-weight: 700; fill: %(GREY)s; }
.ratio { font-size: 12px; font-weight: 700; fill: %(P)s; }
.tick { font-size: 10px; fill: %(GREY)s; }
.setting { left: 28px; top: 318px; width: 262px; font-size: 10.5px; line-height: 1.3; color: %(GREY)s; }
.legend { left: 300px; top: 330px; width: 360px; font-size: 10px; color: %(GREY)s; white-space: nowrap; }
.sw { display: inline-block; width: 11px; height: 11px; vertical-align: -1px; margin: 0 4px 0 0; border-radius: 2px; }
.card { width: 246px; padding: 7px 10px 7px 12px; border-radius: 6px; font-size: 11px; line-height: 1.32; }
.card .t { font-size: 13px; font-weight: 700; margin-bottom: 2px; }
.c1 { left: 668px; top: 80px; background: %(B_PALE)s; border-left: 5px solid %(B)s; }
.c2 { left: 668px; top: 167px; background: %(P_PALE)s; border-left: 5px solid %(P)s; }
.c3 { left: 668px; top: 254px; background: %(O_PALE)s; border-left: 5px solid %(O)s; }
.foot { left: 300px; top: 364px; width: 632px; text-align: right; font-size: 9.5px; color: %(GREY)s; }
</style></head><body>
<div class="abs head">Widely used graviton-exchange rates <span class="em">count the plasma's fermions twice</span></div>
<div class="abs sub">Separate Ways and the Upside Down: corrected rates, and what they mean for the dark matter of a CPT-symmetric universe</div>

<div class="abs" style="left:28px; top:80px;">%(HOUR)s</div>
<div class="abs setting">Boyle, Finn and Turok: the universe before the bang is the CPT image of the one after, and its dark matter is a right-handed neutrino, near 480 PeV, made at the bang.</div>

<div class="abs big">&#215;%(RF)s</div>
<div class="abs bigcap">too large for fermion dark matter, and <span style="color:%(P)s;font-weight:700">&#215;%(RS)s</span> for a scalar: the Standard Model's 45 Weyl fermions counted as 45 Dirac fermions.</div>
<div class="abs" style="left:300px; top:164px;">%(BARS)s</div>
<div class="abs legend"><span class="sw" style="background:%(B)s"></span><b style="color:%(B_DEEP)s">corrected</b>, in Table 1 &nbsp;&nbsp;<span class="sw" style="background:#fff;border:2px solid %(P)s;box-sizing:border-box"></span>published count &nbsp;&nbsp;<span class="tag">Section 2, Fig. 1, Appendix C</span></div>

<div class="abs card c1"><div class="t">Without inflation: %(MC)s PeV at most</div>No starting state that respects CPT makes fewer dark-matter particles than Boyle, Finn and Turok's, so their mass is a ceiling, with no floor. <span class="tag">Section 4</span></div>
<div class="abs card c2"><div class="t">Hot plasma: %(TM)s&#215;10<sup>%(TE)s</sup> GeV is enough</div>Even with the corrected rate, in Einstein gravity a plasma that ever reached that temperature makes all the dark matter. <span class="tag">Section 5</span></div>
<div class="abs card c3"><div class="t">If it decays: a line at %(LINE)s PeV at most</div>A dark-matter neutrino decays to a neutrino line at half its mass. The lightest neutrino is massless, so the masses sum to %(NU)s meV. <span class="tag">Section 6</span></div>

<div class="abs foot">Every number is made by public code and checked against the text: doi.org/10.5281/zenodo.23034297</div>
</body></html>
""" % {"INK": INK, "P": P, "B": B, "GREY": GREY, "B_PALE": B_PALE, "P_PALE": P_PALE, "O": O, "O_PALE": O_PALE, "HOUR": hourglass_svg(),
       "BARS": bars_svg(ratios), "RF": "%.2f" % rf["fermion"], "RS": "%.2f" % rf["scalar"], "MC": "%.1f" % M_ceiling, "B_DEEP": B_DEEP,
       "TM": mant, "TE": int(ex), "LINE": "%.1f" % line, "NU": "%.1f" % nu}

def page_toc(ratios, M_ceiling, line):
    """The same story in three steps and about twenty words, sized to read when ScienceDirect shows the
    image 200 pixels high: every word that carries the message is at least 18 of the page's 384 pixels."""
    rf = dict(ratios)
    SLATE_D, SLATE, SLATE_P = "#334155", "#64748b", "#e2e8f0"
    return """<!doctype html>
<html><head><meta charset="utf-8"><title>Separate Ways v6 graphical abstract, contents-page size</title>
<style>
@page { size: 10in 4in; margin: 0; }
html, body { margin: 0; padding: 0; }
body { width: 960px; height: 384px; position: relative; overflow: hidden; background: #ffffff;
       font-family: Arial, "Liberation Sans", Helvetica, sans-serif; color: %(INK)s; }
.abs { position: absolute; }
.head { left: 28px; top: 14px; width: 904px; font-size: 29px; font-weight: 700; letter-spacing: -0.3px; }
.head .em { color: %(P)s; }
.cap { font-size: 19px; font-weight: 700; line-height: 1.15; text-align: center; }
.big { font-weight: 700; line-height: 1; letter-spacing: -1.5px; text-align: center; }
.side { font-weight: 700; }
.bang { font-size: 13px; font-weight: 700; fill: %(INK)s; }
.cpt { font-size: 13px; font-weight: 700; fill: %(GREY)s; }
.ratio { font-size: 15px; font-weight: 700; fill: %(P)s; }
.tick { font-size: 14px; fill: %(GREY)s; }
.step { font-size: 40px; color: #cbd5e1; font-weight: 700; }
.foot { left: 28px; top: 362px; width: 904px; text-align: right; font-size: 11px; color: %(GREY)s; }
</style></head><body>
<div class="abs head">Widely used graviton-exchange rates <span class="em">count fermions twice</span></div>

<div class="abs" style="left:40px; top:66px;">%(HOUR)s</div>
<div class="abs cap" style="left:20px; top:304px; width:260px;">a CPT-symmetric bang</div>
<div class="abs step" style="left:286px; top:162px;">&#8250;</div>

<div class="abs big" style="left:318px; top:70px; width:330px; font-size:92px; color:%(P)s;">&#215;%(RF)s</div>
<div class="abs" style="left:318px; top:150px;">%(BARS)s</div>
<div class="abs cap" style="left:318px; top:304px; width:330px;">corrected for every kind<br>of gravitational dark matter</div>
<div class="abs step" style="left:660px; top:162px;">&#8250;</div>

<div class="abs big" style="left:690px; top:84px; width:250px; font-size:54px; color:%(B_DEEP)s;">&#8804;%(MC)s<span style="font-size:30px; letter-spacing:0"> PeV</span></div>
<div class="abs cap" style="left:690px; top:146px; width:250px;">dark matter without inflation</div>
<div class="abs big" style="left:690px; top:222px; width:250px; font-size:34px; color:%(O_DEEP)s; letter-spacing:-0.5px;">&#957; line &#8804;%(LINE)s<span style="font-size:22px; letter-spacing:0"> PeV</span></div>
<div class="abs cap" style="left:690px; top:264px; width:250px;">if it decays</div>

<div class="abs foot">Separate Ways and the Upside Down &#183; every number from public code &#183; doi.org/10.5281/zenodo.23034297</div>
</body></html>
""" % {"INK": INK, "P": P, "B": B, "GREY": GREY, "O_DEEP": O_DEEP,
       "HOUR": hourglass_svg(up=(SLATE_P, SLATE_D), down=(SLATE_P, SLATE), side_px=15, width=250),
       "BARS": bars_svg(ratios, toc=True), "RF": "%.2f" % rf["fermion"], "MC": "%.1f" % M_ceiling, "LINE": "%.1f" % line,
       "B_DEEP": B_DEEP}

def render(html_path, stem=STEM):
    base = [CHROME, "--headless=new", "--disable-gpu", "--no-sandbox", "--hide-scrollbars", "--no-pdf-header-footer",
            "--virtual-time-budget=3000"]
    url = "file://" + html_path
    subprocess.run(base + ["--print-to-pdf=" + stem + ".pdf", url], check=True, capture_output=True, timeout=120)
    subprocess.run(base + ["--screenshot=" + stem + ".png", "--window-size=960,384", "--force-device-scale-factor=3.125", url],
                   check=True, capture_output=True, timeout=120)

if __name__ == "__main__":
    ratios, T_star, M_ceiling, line, nu = numbers()
    print("rates: " + ", ".join("%s %.3f" % r for r in ratios) + "; T* = %.3e GeV; ceiling %.1f PeV; line %.1f PeV; sum %.2f meV"
          % (T_star, M_ceiling, line, nu))
    html = STEM + ".html"
    open(html, "w", encoding="utf-8").write(page(ratios, T_star, M_ceiling, line, nu))
    render(html)
    print("wrote graphical_abstract_NPB_v6.html, .pdf and .png")
    html = STEM + "_toc.html"
    open(html, "w", encoding="utf-8").write(page_toc(ratios, M_ceiling, line))
    render(html, STEM + "_toc")
    print("wrote graphical_abstract_NPB_v6_toc.html, .pdf and .png")

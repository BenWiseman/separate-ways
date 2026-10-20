#!/usr/bin/env python3
"""Figures for Separate Ways v6.

fig_v6_occupations.png: occupation n(x) after production, and the integrand x^2 n that fixes the
  mass, for the states of sections 2 and 3 and the state at the upper edge of the family that
  Table A.1 lists at x_c = 0.564. The massless vacuum comes from checks/calc/v6_reheating_production.py.
fig_v6_rates.png: the graviton-exchange coefficients against the published 45-Dirac sums, and Boyle,
  Turok and Vaibhav's Einstein-limit rate corrected one step at a time, read from the printed output
  of checks/calc/v6_grav_rate.py (the same lines v6_numbers_check.py checks), so it adds no number.
"""
import math, os, re, subprocess, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
CALC = os.path.join(HERE, "..", "..", "checks", "calc")
sys.path.insert(0, CALC)
from v6_reheating_production import Profile, n_exact, n_fourier

sys.path.insert(0, HERE)
import v6_style as st
st.apply()

X_C = 1/math.sqrt(math.pi)     # 0.564, the upper-edge state of Table A.1 that gives 330.4 PeV

def occupation(prof, xs):
    return np.array([n_exact(x, prof) if x <= 6.0 else n_fourier(x, prof) for x in xs])

def fig_occupations():
    cache = os.path.join(HERE, "fig_v6_occupations_vac.npz")
    xs = np.concatenate([np.geomspace(0.03, 6.0, 90), np.geomspace(6.0, 300.0, 50)[1:]])
    if os.path.exists(cache):
        n_vac = np.load(cache)["n_vac"]
    else:
        n_vac = occupation(Profile(t_rh=1e-4), xs)      # the massless spinor at t = 0, projected onto the out states
        np.savez(cache, n_vac=n_vac)
    # the closed-form states on a finer grid that has the jump at x_c on it
    xf = np.sort(np.concatenate([np.geomspace(0.03, 300.0, 1200), [X_C*(1 - 1e-9), X_C*(1 + 1e-9)]]))
    n_half = (1 - np.sqrt(-np.expm1(-math.pi*xf**2)))/2
    n_in = np.exp(-math.pi*xf**2)
    n_edge = np.where(xf < X_C, 1 - n_half, n_half)
    fig, (a, b) = plt.subplots(1, 2, figsize=(st.TEXT_WIDTH, 2.8))
    curves = [(xf, n_in, st.GREY, ":", 1.3, "in-vacuum, not invariant"),
              (xf, n_half, st.BLUE, "-", 1.9, "half-angle state (491.6 PeV)"),
              (xf, n_edge, st.PURPLE, "-.", 1.4, r"upper edge below $x_c$ (330.4 PeV)"),
              (xs, n_vac, st.INK, "--", 1.1, "massless vacuum at the bang")]
    for x, n, c, ls, lw, lab in curves:
        ok = n > 1e-14
        z = 3 if c == st.BLUE else 2                     # the ceiling's state on top where the edge state joins it
        a.loglog(x[ok], n[ok], color=c, ls=ls, lw=lw, label=lab, zorder=z)
        sel = x <= 3.0
        b.plot(x[sel], (x**2*n)[sel], color=c, ls=ls, lw=lw, zorder=z)
    a.loglog(xs[xs > 3], 1/(16*xs[xs > 3]**4), color=st.GREY_L, lw=0.8, ls="-", zorder=0)
    a.text(30, 1/(16*30**4)*6, r"$1/16x^4$", color=st.GREY, fontsize=8.5)
    a.set_xlabel(r"$x=p/\sqrt{Ma_1}$"); a.set_ylabel(r"occupation $n(x)$")
    a.set_ylim(1e-9, 2); a.set_xlim(0.03, 300)
    b.set_xlabel(r"$x$"); b.set_ylabel(r"$x^2 n(x)$")
    b.set_xlim(0, 3); b.set_ylim(0, 0.34)
    st.panel(a, "(a)"); st.panel(b, "(b)")
    h, l = a.get_legend_handles_labels()
    order = [1, 2, 0, 3]                                 # the ceiling first, then the state that undercuts it
    b.legend([h[i] for i in order], [l[i] for i in order], loc="upper right", fontsize=7.4, handlelength=2.4,
             borderaxespad=0.3, labelspacing=0.35)
    fig.tight_layout(w_pad=1.6)
    fig.savefig(os.path.join(HERE, "fig_v6_occupations.png"))
    fig.savefig(os.path.join(HERE, "fig_v6_occupations.pdf"))
    plt.close(fig)
    # the plotted upper-edge state must be Table A.1's: I/I_* = 1 + int_0^x_c x^2 (1 - 2 n_*) dx / (sqrt(pi) I_*)
    x = np.linspace(0, X_C, 200001); nh = (1 - np.sqrt(1 - np.exp(-math.pi*x**2)))/2
    I_star = 0.0127596674
    return 1 + np.trapezoid(x**2*(1 - 2*nh), x)/math.sqrt(math.pi)/I_star

def fig_rates():
    """(a) corrected coefficients against the 45-Dirac sums; (b) BTV's Einstein limit, step by step.
    Every value is parsed from v6_grav_rate.py's printed output, the lines the numbers check reads."""
    out = subprocess.run([sys.executable, os.path.join(CALC, "v6_grav_rate.py")], capture_output=True,
                         text=True, timeout=1200).stdout
    rows = []
    for lab, key in (("scalar,\n$\\xi=0$", r"real scalar, xi = 0"), ("scalar,\n$\\xi=1/6$", r"real scalar, xi = 1/6"),
                     ("Dirac\nfermion", r"Dirac fermion"), ("massless\nvector", r"massless vector")):
        m = re.search(re.escape(key) + r"\s+r = .*?R = ([0-9.e+-]+);.*?with 45 Dirac ([0-9.e+-]+), ([0-9.]+) times", out)
        if not m: raise SystemExit("v6_grav_rate.py output lacks the %s row" % key)
        rows.append((lab.replace("\\n", "\n"), float(m.group(1)), float(m.group(2)), float(m.group(3))))
    alpha = rows[2][1]
    m1 = re.search(r"BTV Einstein limit as printed: ([0-9.e+-]+)", out)
    m2 = re.search(r"transverse gauge vertex: ([0-9.e+-]+); Majorana pair counting: ([0-9.e+-]+); identical real scalars: ([0-9.e+-]+)", out)
    m3 = re.search(r"with 45 Weyl fermions in place of their 47: ([0-9.e+-]+)", out)
    if not (m1 and m2 and m3): raise SystemExit("v6_grav_rate.py output lacks the BTV steps")
    steps = [("as\nprinted", float(m1.group(1))), ("gauge\nvertex", float(m2.group(1))),
             ("one $N_1$\nper event", float(m2.group(2))), ("1/2 for\nscalars", float(m2.group(3))),
             ("45 Weyl,\nnot 47", float(m3.group(1))), ("this\npaper", alpha)]
    fig, (a, b) = plt.subplots(1, 2, figsize=(st.TEXT_WIDTH, 3.0), gridspec_kw={"width_ratios": [1, 1.12]})
    xs = np.arange(len(rows)); w = 0.36
    a.bar(xs - w/2, [r[1] for r in rows], w, color=st.BLUE, label="corrected (Table 1)")
    a.bar(xs + w/2, [r[2] for r in rows], w, color="white", edgecolor=st.PURPLE, lw=1.3,
          label="published count (45 Dirac fermions)")
    for x, r in zip(xs, rows):
        a.text(x + w/2, r[2]*1.13, r"$\times%.2f$" % r[3], ha="center", va="bottom", fontsize=8.2, color=st.PURPLE)
    a.set_yscale("log"); a.set_ylim(1e-3, 1.5e-1)
    a.set_xticks(xs); a.set_xticklabels([r[0] for r in rows], fontsize=8.2)
    a.tick_params(axis="x", which="both", top=False, bottom=False)
    a.set_ylabel(r"rate from the plasma ($T^8/M_P^4$)")
    a.legend(loc="upper left", fontsize=7.8)
    ys = [v/alpha for _, v in steps]
    b.step(range(len(steps)), ys, where="mid", color=st.GREY, lw=1.0)
    b.plot(range(len(steps) - 1), ys[:-1], "o", color=st.PURPLE, ms=5.5, zorder=3)
    b.plot([len(steps) - 1], [1.0], "o", color=st.BLUE, ms=6.5, zorder=3)
    for i, y in enumerate(ys):
        lab = "1" if abs(y - 1) < 1e-9 else ("%.2f" % y if y >= 2 else "%.3f" % y)   # as Table B.1 prints them
        b.text(i, y + 0.22, lab, ha="center", va="bottom", fontsize=8.2)
    b.set_xticks(range(len(steps))); b.set_xticklabels([s_ for s_, _ in steps], fontsize=7.2)
    b.tick_params(axis="x", which="both", top=False, bottom=False)
    b.set_xlim(-0.5, len(steps) - 0.5); b.set_ylim(0, 5.2)
    b.set_ylabel("Boyle, Turok and Vaibhav's rate\n" + r"in units of $\alpha$")
    st.panel(a, "(a)"); st.panel(b, "(b)")
    fig.tight_layout(w_pad=1.4)
    fig.savefig(os.path.join(HERE, "fig_v6_rates.png"))
    fig.savefig(os.path.join(HERE, "fig_v6_rates.pdf"))
    plt.close(fig)
    return rows, ys

if __name__ == "__main__":
    r = fig_occupations()
    print("occupations figure: the upper-edge state below x_c = %.3f has I/I_* = %.3f, Table A.1's 2.700: %s"
          % (X_C, r, "yes" if abs(r - 2.700) < 1e-3 else "NO  <-- FAILED"))
    if abs(r - 2.700) >= 1e-3: sys.exit(1)
    rows, ys = fig_rates()
    print("rates figure: " + "; ".join("%s %.4e vs %.4e (x%.3f)" % (r[0].replace(chr(10), " "), r[1], r[2], r[3]) for r in rows))
    print("rates figure: BTV steps in units of alpha " + ", ".join("%.3f" % y for y in ys))

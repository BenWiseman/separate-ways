#!/usr/bin/env python3
"""Figures for "Upside Down Inflation".

fig_infl_states.png: occupation n(x) and the integrand x^2 n after inflation at H_inf = 1e12 GeV, with a
  sudden end of reheating at two temperatures, the massless vacuum (instant reheating), the half-angle
  state for comparison, and the bound over CPT-invariant states averaged over the interference phase.
  Computed with checks/calc/v6_reheating_production.py.
fig_infl_mass.png: the N1 mass against the reheating temperature, from the mass tables that
  v6_mass_with_bath.py and v6_mass_smooth.py write in checks/calc.
fig_inflaton.png: (a) the T-model potential, where the bang may sit, where today's Hubble scale leaves
  and where inflation ends, with the factor dchi/dphi that shields the plateau; (b) the N1 mass against
  the reheating temperature, with the reach of the two couplings. Its numbers are parsed from the
  printed output of checks/calc/v6_alpha_inflaton.py.
"""
import math, os, re, subprocess, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
CALC = os.path.join(HERE, "..", "..", "checks", "calc")
sys.path[:0] = [HERE, os.path.join(HERE, "..", "paper2"), CALC]   # in the release v6_style.py sits beside this file
from v6_reheating_production import Profile, n_exact, n_fourier, n_star_avg
import v6_style as st
st.apply()

MP, GSTAR = 2.435e18, 106.75
CH = math.sqrt(math.pi**2*GSTAR/90)

def profile_for(T_rh, M, H_inf=1e12):
    H_rh = CH*T_rh**2/MP
    t_rh = math.sqrt(M/H_rh); t_e = min((8*t_rh*M/H_inf)**(1/3) - t_rh, t_rh)
    return Profile(t_rh=t_rh, t_e=t_e)

def occupation(prof, xs):
    return np.array([n_exact(x, prof) if x <= 6.0 else n_fourier(x, prof) for x in xs])

def fig_states():
    cache = os.path.join(HERE, "fig_infl_states.npz")
    xs = np.concatenate([np.geomspace(0.03, 6.0, 90), np.geomspace(6.0, 300.0, 50)[1:]])
    if os.path.exists(cache):
        d = np.load(cache); n_sud, n_18, n_50 = d["n_sud"], d["n_18"], d["n_50"]
    else:
        n_sud = occupation(Profile(t_rh=1e-4), xs)
        n_18 = occupation(profile_for(1.87e13, 4.134e8), xs)
        n_50 = occupation(profile_for(5e13, 2.784e8), xs)
        np.savez(cache, n_sud=n_sud, n_18=n_18, n_50=n_50)
    n_half = (1 - np.sqrt(1 - np.exp(-math.pi*xs**2)))/2
    fig, (a, b) = plt.subplots(1, 2, figsize=(st.TEXT_WIDTH, 2.9))
    curves = [(n_half, st.PURPLE, "-", 1.4, "half-angle state, no inflation"),
              (n_sud, st.INK, "--", 1.1, "massless vacuum, instant reheating"),
              (n_18, st.BLUE, "-", 1.7, r"$T_{\rm rh}=1.87\times10^{13}$ GeV, sudden end"),
              (n_50, st.BLUE_M, "-.", 1.3, r"$T_{\rm rh}=5\times10^{13}$ GeV, sudden end"),
              (n_star_avg(n_18), st.ORANGE, "--", 1.4, r"bound over CPT-invariant states, $1.87\times10^{13}$ GeV")]
    for n, c, ls, lw, lab in curves:
        ok = n > 1e-14
        a.loglog(xs[ok], n[ok], color=c, ls=ls, lw=lw, label=lab)
        sel = xs <= 3.0
        b.plot(xs[sel], (xs**2*n)[sel], color=c, ls=ls, lw=lw)
    a.set_xlabel(r"$x=p/\sqrt{Ma_1}$"); a.set_ylabel(r"occupation $n(x)$")
    a.set_ylim(1e-9, 2); a.set_xlim(0.03, 300)
    b.set_xlabel(r"$x$"); b.set_ylabel(r"$x^2 n(x)$")
    b.set_xlim(0, 3); b.set_ylim(0, 0.13)
    st.panel(a, "(a)"); st.panel(b, "(b)")
    h, l = a.get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", ncol=2, bbox_to_anchor=(0.5, -0.02), columnspacing=1.6)
    fig.tight_layout(rect=(0, 0.2, 1, 1), w_pad=1.6)
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(HERE, "fig_infl_states.%s" % ext))
    plt.close(fig)

def read_tsv(name):
    rows = {}
    with open(os.path.join(CALC, name)) as f:
        next(f)
        for line in f:
            parts = line.rstrip("\n").split("\t")
            tag, T, H, M, fb = parts[:5]
            rows.setdefault(tag, []).append((float(T), float(M), float(fb)))
    return {k: np.array(sorted(v)) for k, v in rows.items()}

def fig_mass():
    d = read_tsv("v6_mass_with_bath.tsv")             # sudden end of reheating
    p = read_tsv("v6_mass_smooth.tsv")                # perturbative reheating, quadratic minimum
    q = read_tsv("v6_mass_smooth_w13.tsv")            # perturbative reheating, quartic minimum
    pec = read_tsv("v6_mass_smooth_ec.tsv")["with plasma"]       # the same with the minimal Einstein-Cartan term
    qec = read_tsv("v6_mass_smooth_w13_ec.tsv")["with plasma"]
    fig, ax = plt.subplots(figsize=(st.TEXT_WIDTH*0.82, 3.9))
    su, pe, pw, pb, qw = (d["central: alpha 1.2e-2, no boost"], p["expansion only"], p["with plasma"],
                          p["bound, with plasma"], q["with plasma"])
    Tg = np.geomspace(pw[0, 0], pw[-1, 0], 200)
    lo = np.exp(np.interp(np.log(Tg), np.log(qw[:, 0]), np.log(qw[:, 1])))
    hi = np.exp(np.interp(np.log(Tg), np.log(pw[:, 0]), np.log(pw[:, 1])))
    ax.fill_between(Tg, np.minimum(lo, hi)/1e6, np.maximum(lo, hi)/1e6, color=st.SHADE, lw=0,
                    label="between quadratic and quartic minima")
    ax.plot(pw[:, 0], pw[:, 1]/1e6, color=st.BLUE, lw=1.9, label="quadratic minimum, with graviton exchange")
    ax.plot(pe[:, 0], pe[:, 1]/1e6, color=st.GREY, ls="--", lw=1.2, label="quadratic minimum, expansion alone")
    ax.plot(su[:, 0], su[:, 1]/1e6, color=st.BLUE_M, lw=1.1, ls="-.", label="quadratic, sudden end of reheating")
    ax.plot(pec[:, 0], pec[:, 1]/1e6, color=st.BLUE, lw=1.2, ls=":", label="quadratic, with the Einstein-Cartan term")
    ax.plot(qw[:, 0], qw[:, 1]/1e6, color=st.PURPLE, lw=1.9, label="quartic minimum, with graviton exchange")
    ax.plot(qec[:, 0], qec[:, 1]/1e6, color=st.PURPLE, lw=1.2, ls=":", label="quartic, with the Einstein-Cartan term")
    # where H = M falls at T_rh on the quadratic perturbative curve
    Tb = math.exp(np.interp(0.0, np.log(pw[:, 1]/(CH*pw[:, 0]**2/MP))[::-1], np.log(pw[:, 0])[::-1]))
    ax.axvline(Tb, color=st.GREY_L, lw=0.9, zorder=0)
    ax.text(Tb*0.93, 1450, "made during\nreheating", ha="right", va="center", fontsize=8.2, color=st.GREY)
    ax.text(Tb*1.07, 1450, "made in the\nradiation era", ha="left", va="center", fontsize=8.2, color=st.GREY)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlim(2.5e12, 2.2e14); ax.set_ylim(5, 2500)       # the Einstein-Cartan curves reach 6.6 PeV
    ax.set_xlabel(r"reheating temperature $T_{\rm rh}$ (GeV)")
    ax.set_ylabel(r"$N_1$ mass $M$ (PeV)")
    h, l = ax.get_legend_handles_labels()
    order = [1, 2, 3, 4, 5, 6, 0]
    fig.legend([h[i] for i in order], [l[i] for i in order], loc="lower center", ncol=2, bbox_to_anchor=(0.5, -0.01),
               columnspacing=1.4)
    fig.tight_layout(rect=(0, 0.2, 1, 1))
    fig.savefig(os.path.join(HERE, "fig_infl_mass.png"))
    fig.savefig(os.path.join(HERE, "fig_infl_mass.pdf"))
    plt.close(fig)
    return Tb

def grab(pattern, text, what):
    m = re.search(pattern, text)
    if not m: raise SystemExit("v6_alpha_inflaton.py output lacks %s" % what)
    return m

def fig_model():
    r = subprocess.run([sys.executable, os.path.join(CALC, "v6_alpha_inflaton.py")], capture_output=True, text=True, timeout=3600)
    if r.returncode != 0: raise SystemExit("v6_alpha_inflaton.py exited %d" % r.returncode)
    out = r.stdout
    alpha = float(grab(r"alpha = ([0-9.]+), r =", out, "alpha").group(1))
    phi_end = float(grab(r"inflation ends \(eps = 1\) at phi = ([0-9.]+) M_P", out, "the end of inflation").group(1))
    phi58 = float(grab(r"leaves at phi = ([0-9.]+) M_P", out, "phi(N = 58)").group(1))
    T_tenth = float(grab(r"c is a tenth of V0 at y = [0-9.]+, T_rh = ([0-9.e+]+) GeV", out, "the tenth").group(1))
    T_lim = float(grab(r"and equals V0 at y = [0-9.]+, T_rh = ([0-9.e+]+) GeV", out, "the limit").group(1))
    M_lim = float(grab(r"c at V0: T_rh = [0-9.e+]+ GeV, M = ([0-9]+) PeV", out, "the mass at the limit").group(1))
    T_inst = float(grab(r"above instant reheating, ([0-9.e+]+) GeV", out, "instant reheating").group(1))

    fig, (a, b) = plt.subplots(1, 2, figsize=(st.TEXT_WIDTH, 2.9))
    s = math.sqrt(6*alpha)
    ph = np.linspace(0, 1.3, 600)
    a.plot(ph, np.tanh(ph/s)**2, color=st.INK, lw=1.6)
    a.text(0.36, 1.03, r"$V/V_0=\tanh^2(\phi/\sqrt{6\alpha}\,M_P)$", fontsize=8.5, va="bottom")
    a.axvspan(phi58, 1.3, color=st.BLUE_P, lw=0, zorder=0)
    a.text(1.08, 0.55, "the bang,\nfield at rest", ha="center", fontsize=8.2, color=st.BLUE_D)
    for x, lab in ((phi_end, "inflation\nends"), (phi58, "today's Hubble\nscale leaves")):
        a.axvline(x, color=st.GREY_L, lw=0.9, zorder=0)
        a.text(x - 0.02, 0.25, lab, ha="right", fontsize=8.2, color=st.GREY)
    a.set_xlabel(r"$\phi/M_P$"); a.set_ylabel(r"$V/V_0$"); a.set_xlim(0, 1.3); a.set_ylim(0, 1.15)
    a2 = a.twinx()
    a2.semilogy(ph, 1/np.cosh(ph/s)**2, color=st.PURPLE, lw=1.3, ls="-.")
    a2.set_ylabel(r"$d\chi/d\phi$", color=st.PURPLE); a2.set_ylim(1e-7, 2)
    a2.tick_params(axis="y", which="both", colors=st.PURPLE)
    a.tick_params(axis="y", which="both", right=False)
    st.panel(a, "(a)")

    rows = [l.rstrip("\n").split("\t") for l in open(os.path.join(CALC, "v6_mass_smooth.tsv")) if l.startswith("with plasma\t")]
    T = np.array([float(x[1]) for x in rows]); M = np.array([float(x[3]) for x in rows])/1e6; o = np.argsort(T); T, M = T[o], M[o]
    b.loglog(T, M, color=st.BLUE, lw=1.8)
    b.axvspan(T[0], T_lim, color=st.BLUE_P, lw=0, zorder=0)
    b.axvspan(T_lim, T[-1], color=st.ORANGE_P, lw=0, zorder=0)
    b.axvline(T_tenth, color=st.BLUE_D, lw=0.9, ls=":")
    b.plot([T_lim], [M_lim], "o", color=st.BLUE_D, ms=5, zorder=3)
    b.text(T_lim*0.93, 25, "Yukawa through $\\chi$", ha="right", fontsize=8.2, color=st.BLUE_D)
    b.text(T_lim*1.12, 24, "trilinear through $\\chi$:\nto instant reheating,\nwith tachyonic preheating",
           ha="left", va="bottom", fontsize=7.8, color=st.ORANGE_D)
    b.text(T_lim*1.05, M_lim*1.05, "%d PeV" % M_lim, fontsize=8.2, color=st.BLUE_D, va="bottom")
    b.set_xlabel(r"reheating temperature $T_{\rm rh}$ (GeV)"); b.set_ylabel(r"$N_1$ mass (PeV)")
    b.set_xlim(T[0], T[-1]); b.set_ylim(20, 2000)
    b.text(T[-1]*0.85, 1500, "quadratic minimum,\nwith graviton exchange", ha="right", va="top", fontsize=8.2, color=st.BLUE_D)
    st.panel(b, "(b)")
    fig.tight_layout(w_pad=1.6)
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(HERE, "fig_inflaton.%s" % ext))
    return ("alpha %.4f, phi_end %.2f, phi58 %.2f, T_rh limits %.2e and %.2e GeV, M %d PeV at the limit, instant %.1e GeV"
          % (alpha, phi_end, phi58, T_tenth, T_lim, M_lim, T_inst))

if __name__ == "__main__":
    fig_states()
    print(fig_model())
    Tb = fig_mass()
    print("figures written; H(T_rh) = M on the perturbative curve at T_rh = %.3e GeV" % Tb)
    for name, tag in (("v6_mass_smooth.tsv", "with plasma"), ("v6_mass_with_bath.tsv", "central: alpha 1.2e-2, no boost"),
                      ("v6_mass_smooth_w13.tsv", "with plasma")):
        ce = read_tsv(name)[tag]
        lt, lm = np.log(ce[:, 0]), np.log(ce[:, 1])
        if lm.max() < math.log(4.4e8):
            print("%s: never reaches 440 PeV (largest %.1f PeV)" % (name, math.exp(lm.max())/1e6)); continue
        T440 = math.exp(np.interp(math.log(4.4e8), lm[::-1], lt[::-1]))
        print("%s: M = 440 PeV (line 220 PeV) at T_rh = %.2e GeV" % (name, T440))
    pw = read_tsv("v6_mass_smooth.tsv")["with plasma"]; pb = read_tsv("v6_mass_smooth.tsv")["bound, with plasma"]
    su = read_tsv("v6_mass_with_bath.tsv")["central: alpha 1.2e-2, no boost"]
    lt, lm, fs = np.log(pw[:, 0]), np.log(pw[:, 1]), pw[:, 2]
    print("perturbative, with plasma: under 300 PeV above T_rh = %.2e GeV" % math.exp(np.interp(math.log(3e8), lm[::-1], lt[::-1])))
    print("perturbative, with plasma: plasma makes half the dark matter at T_rh = %.2e GeV" % math.exp(np.interp(0.5, fs, lt)))
    print("bound over CPT-invariant states / massless vacuum (perturbative, with plasma): up to %.1f per cent higher"
          % (100*(max(pb[:, 1]/pw[:, 1]) - 1)))
    i = int(np.argmin(abs(pw[:, 0]/1.87e13 - 1))); j = int(np.argmin(abs(su[:, 0]/1.87e13 - 1)))
    print("at 1.87e13 GeV a sudden end of reheating gives %.1f per cent less than perturbative reheating" % (100*(1 - su[j, 1]/pw[i, 1])))
    qw = read_tsv("v6_mass_smooth_w13.tsv")["with plasma"]; qe = read_tsv("v6_mass_smooth_w13.tsv")["expansion only"]
    k = int(np.argmin(abs(qw[:, 0]/1.87e13 - 1)))
    print("quartic minimum: expansion alone %.1f to %.1f PeV over the grid; with plasma %.1f PeV at 1.87e13 GeV"
          % (qe[:, 1].min()/1e6, qe[:, 1].max()/1e6, qw[k, 1]/1e6))
    lq, mq = np.log(qw[:, 0]), np.log(qw[:, 1])
    print("quartic minimum, with plasma: under 250 PeV above T_rh = %.2e GeV; plasma makes half at %.2e GeV"
          % (math.exp(np.interp(math.log(2.5e8), mq[::-1], lq[::-1])), math.exp(np.interp(0.5, qw[:, 2], lq))))
    for name, lab in (("v6_mass_smooth_ec.tsv", "quadratic"), ("v6_mass_smooth_w13_ec.tsv", "quartic")):
        e = read_tsv(name)["with plasma"]; le = np.log(e[:, 0])
        m = int(np.argmin(abs(e[:, 0]/1.87e13 - 1)))
        print("Einstein-Cartan, %s minimum: %.1f PeV at 1.87e13 GeV (plasma %.1f%%); plasma makes half at %.2e GeV"
              % (lab, e[m, 1]/1e6, 100*e[m, 2], math.exp(np.interp(0.5, e[:, 2], le))))
    # how far apart the sudden and perturbative treatments of a quadratic minimum put T_rh for one mass
    lp, mp_ = np.log(pw[:, 0]), np.log(pw[:, 1]); ls_, ms_ = np.log(su[:, 0]), np.log(su[:, 1])
    worst = 0.0
    for Mx in np.geomspace(1.0e8, 1.0e9, 60):
        if not (mp_.min() < math.log(Mx) < mp_.max() and ms_.min() < math.log(Mx) < ms_.max()): continue
        tp = math.exp(np.interp(math.log(Mx), mp_[::-1], lp[::-1])); ts = math.exp(np.interp(math.log(Mx), ms_[::-1], ls_[::-1]))
        worst = max(worst, tp/ts - 1)
    print("quadratic minimum: sudden and perturbative reheating read one mass (100 to 1000 PeV) as T_rh up to %.0f per cent apart" % (100*worst))

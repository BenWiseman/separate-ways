#!/usr/bin/env python3
"""Every number "Upside Down Inflation" quotes must be in the output of the script that makes it.

    python3 pub/inflaton/inflaton_check.py              # check
    python3 pub/inflaton/inflaton_check.py --selftest   # planted wrong numbers must fail

Each claim is (text in the paper, source, text in that source's output). A source is a script, run as
Separate Ways v6's own check runs it (v6_numbers_check.run), from this directory if it sits here and
from the calc directory otherwise, or "log:" and a file of the unpublished clock test, which exists in
the working tree but not in the code release; there those claims are reported as skipped, by name, and
never silently passed. Table 1 is checked cell by cell against mass_table.py, which reads the solvers'
mass tables. Literature values (Planck, ACT, Boyle, Finn and Turok's 480 PeV) are checked against their
sources elsewhere and are not listed here.
"""
import io, os, re, sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path[:0] = [HERE, os.path.join(HERE, "..", "paper2")]     # in the release v6_numbers_check.py sits beside this file
from v6_numbers_check import run, flat, ROOT, CALC, TSVS, REC, current, sources

NOTE = os.path.join(HERE, "INFLATON_v1.md")
CLOCK = os.path.join(ROOT, "beyond", "clock_expansion", "logs")
A, T, B = "v6_alpha_inflaton.py", "v6_trilinear_bath.py", "v6_bath.py"

CLAIMS = [
 # abstract and Section 3
 (r"A T-model $\alpha$-attractor with $\alpha=0.0045$", A, "alpha = 0.0045, r = 1.60e-05"),
 (r"with $r=1.6\times10^{-5}$ and $n_s=0.9655$", A, "r = 1.60e-05, n_s = 0.9655"),
 (r"which moves it by $2.0\times10^{-5}M_P$", A, "moves the field by -2.0e-05 M_P"),
 (r"reheats it to at most $1.2\times10^{13}$ GeV, which puts the dark matter at 561 PeV or heavier", A, "c at V0: T_rh = 1.18e+13 GeV, M = 561 PeV"),
 (r"$m_\phi=H_{\rm inf}/\sqrt\alpha=1.5\times10^{13}$ GeV", A, "inflaton mass at the minimum 1.49e+13 GeV"),
 (r"Today's Hubble scale leaves the horizon at $\phi=0.86\,M_P$, and inflation ends at $0.29\,M_P$", A, "today's Hubble scale leaves at phi = 0.86 M_P"),
 (r"Today's Hubble scale leaves the horizon at $\phi=0.86\,M_P$, and inflation ends at $0.29\,M_P$", A, "inflation ends (eps = 1) at phi = 0.29 M_P"),
 (r"lies $2.8\sigma$ below the $0.974\pm0.003$", A, "and 2.8 sigma below ACT DR6's 0.974 +- 0.003"),
 (r"moves the field only $1.6\times10^{-4}M_P$ from $0.86\,M_P$, 0.11 e-folds", A, "the field's own slope moves it -1.6e-04 M_P, -0.11 e-folds"),
 (r"at a temperature of $8.4\times10^{14}$ GeV", T, "the plasma meets the plateau at T = 8.44e+14 GeV"),
 # Section 4
 (r"$d\chi/d\phi=\mathrm{sech}^2u$ is $1.2\times10^{-4}$ at $0.86\,M_P$", A, "sech^2 u there = 1.2e-04"),
 (r"$\sqrt{6\alpha}\,M_P=0.16\,M_P$", A, "sqrt(6 alpha) M_P = 0.16 M_P"),
 (r"the plasma moves the field $2.0\times10^{-5}M_P$, which shortens inflation by 0.014 e-folds", A, "moves the field by -2.0e-05 M_P, which changes the e-folds of inflation by -1.4e-02"),
 (r"never exceeds $1.7\times10^{-4}H^2$", A, "largest |dF/dphi|/H^2 1.7e-04"),
 (r"drags the field from 0.86 to $0.08\,M_P$", A, "the field ends at 0.08 M_P from 0.86"),
 # Section 5
 (r"$T_{\rm rh}=y\times6.5\times10^{14}$ GeV", A, "T_rh = y x 6.50e+14 GeV"),
 (r"holding $c$ to a tenth of $V_0$ gives $y=0.0102$ and $T_{\rm rh}=6.7\times10^{12}$ GeV", A, "c is a tenth of V0 at y = 0.0102, T_rh = 6.65e+12 GeV"),
 (r"the potential rises all the way to its plateau only while $c<V_0$", A, "plant: the potential rises to the plateau for c = 0.99 V0 and not for 1.01 V0      yes"),
 (r"That allows $y$ up to 0.0182, where $T_{\rm rh}=1.2\times10^{13}$ GeV", A, "and equals V0 at y = 0.0182, T_rh = 1.18e+13 GeV"),
 (r"Instant reheating, at $8.4\times10^{14}$ GeV, needs $\mu=1.9\times10^{13}$ GeV", A, "instant reheating needs mu = 1.9e+13 GeV"),
 (r"Instant reheating, at $8.4\times10^{14}$ GeV, needs $\mu=1.9\times10^{13}$ GeV", A, "above instant reheating, 8.4e+14 GeV"),
 (r"is 3.2 per cent of the plateau's and its slope 2.1 per cent of the tree-level slope", A, "is 3.2 per cent of the plateau's and its slope 2.1 per cent of the tree-level slope"),
 (r"moves the field $7.2\times10^{-6}M_P$", A, "moves the field by -7.2e-06 M_P"),
 # Section 6
 (r"gives 1112 PeV at $3\times10^{12}$ GeV, 681 PeV at $8\times10^{12}$ GeV", "mass_table.py", r"| $3\times10^{12}$ | 1112 | 1112 |"),
 (r"gives 1112 PeV at $3\times10^{12}$ GeV, 681 PeV at $8\times10^{12}$ GeV", "mass_table.py", r"| $8\times10^{12}$ | 681 | 681 |"),
 (r"or 747 PeV with $c$ held to a tenth of", A, "c at a tenth of V0: T_rh = 6.65e+12 GeV, M = 747 PeV"),
 (r"above the 491.6 PeV that bounds it without inflation", "v6_abundance.py", "M_* = 491.6 PeV"),
 # Section 7
 (r"$2/3\pi\approx0.21$ new records per Hubble volume per Hubble time", A, "2/(3 pi) = 0.21 per Hubble volume per Hubble time"),
 (r"by 9.8 per cent", "log:clock_de_h0.log", "shift +9.83 per cent"),
 (r"by 83.5 in $\chi^2$", "log:clock_de_h0.log", "chi2 - chi2(LCDM) = +83.49: FAILS"),
 (r"a free amplitude does not rescue it", "log:clock_de_followups.log", "B1(f = 0.1): H0 = 68.62, omega_m = 0.1414, dchi2 vs LCDM +0.47"),
 # the model's own checks
 (r"Section 4 shows that couplings through the inflaton itself fail", A, "all checks pass"),
 # Sections 2 to 4, Appendix A and the conclusions, from the sections Separate Ways v6 had on inflation
 ("Yukawa coupling to $N_1$ is below about $M/20M_P$", "v6_bath.py", "so y below 0.05 M/M_P keeps it under a tenth"),
 ("$\\epsilon=1.0\\times10^{-6}$ and gives $r=1.6\\times10^{-5}$", "v6_bath.py", "r = 1.6e-05, slow-roll epsilon = 1.0e-06"),
 ("left the horizon 58 e-folds", "v6_bath.py", "left the horizon 58.2 e-folds"),
 ("needs $y=0.035$", "v6_bath.py", "needs y = 0.035"),
 ("mass of $1.2\\times10^{13}$ GeV, 12 times $H_{\\rm inf}$, when the plasma at $8.4\\times10^{14}$ GeV", "v6_bath.py", "at T = 8.4e+14 GeV its thermal mass is 1.2e+13 GeV, 12 H_inf"),
 ("with yield $2.5\\times10^{-5}$", "v6_bath.py", "2.50e-05 (x2.8e+13"),
 ("so 12.6 e-folds bring it", "v6_bath.py", "N > 12.6"),
 ("exceeds $10^6$ already at $N=12.6$", "v6_state_at_bang.py", "x at the bang = 1.06e+06"),
 ("$5\\times10^{-26}$ at $x=1.06\\times10^6$", "v6_state_at_bang.py", "1/(16 x^4) = 5.0e-26"),
 ("about $10^{22}$ for 60 e-folds", "v6_state_at_bang.py", "k/(aT) = 8.0e+21"),
 ("once inflation lasts 12.7 to 14.1 e-folds", "v6_state_at_bang.py", "N = 12.7 (bound), 13.2 (up), 13.8 (down)"),
 ("once inflation lasts 12.7 to 14.1 e-folds", "v6_state_at_bang.py", "N = 12.9 (bound), 13.4 (up), 14.1 (down)"),
 ("at 58 e-folds $10^{78}$ times as much", "v6_bath.py", "carries 1e+78 times the plasma's energy"),
 ("to $7.7\\times10^{13}$ GeV when", "v6_mass_smooth.py", "T_max = 7.67e+13 GeV"),
 ("it is 261.8 PeV at every", "v6_mass_smooth.py --w13", "expansion-only mass from 261.8 to 261.8 PeV"),
 ("269.4, 261.3 and 258.8 PeV for", "v6_instant.py", "M = 269.4 PeV"),
 ("269.4, 261.3 and 258.8 PeV for", "v6_instant.py", "M = 261.3 PeV"),
 ("269.4, 261.3 and 258.8 PeV for", "v6_instant.py", "M = 258.8 PeV"),
 ("449 PeV at $T_{\\rm rh}=1.87\\times10^{13}$ GeV, 1112", "v6_mass_smooth.py", "1.87e+13       449.4"),
 ("449 PeV at $T_{\\rm rh}=1.87\\times10^{13}$ GeV, 1112", "v6_mass_smooth.py", "3.00e+12      1111.6      555.8        0.0%"),
 ("makes 0.6 per cent of the abundance", "v6_mass_smooth.py", "449.4      224.7        0.6%"),
 ("8 per cent at $5\\times10^{13}$ GeV", "v6_mass_smooth.py", "304.2      152.1        8.0%"),
 ("most of it above $1.05\\times10^{14}$ GeV", "fig_inflaton.py", "half the dark matter at T_rh = 1.05e+14 GeV"),
 ("runs from $4.2\\times10^5$ GeV", "v6_bath.py", "the plasma sets M = 4.2e+05 GeV"),
 ("near $3.7\\times10^6$ GeV", "v6_bath.py", "at T_rh = 3.7e+06 GeV"),
 ("differ by 3 to 10 per cent in the mass", "v6_reheating_smooth.py", "1.031"),
 ("differ by 3 to 10 per cent in the mass", "v6_reheating_smooth.py", "1.102"),
 ("is 0.53 of the sudden formula", "v6_mass_smooth.py", "ratio 0.53; T_max = 7.67e+13"),
 ("0.82 of ours for the same", "v6_bath.py", "Garcia et al. T_RH / ours = 0.816"),
 ("by at most 0.5 per cent between", "v6_hinf.py", "0.50% for H_inf = 1e11 or 1e13"),
 ("it is 3 per cent higher at", "v6_hinf.py", "1.87e+13        424.8"),
 ("10 per cent higher at $3\\times10^{12}$ GeV", "v6_hinf.py", "10.1% for H_inf = 1e10"),
 ("from about $9\\times10^{13}$ GeV for a constant rate to $3.1\\times10^{13}$ and $1.1\\times10^{14}$ GeV", "v6_quartic_decay.py", "plasma makes half at 9.11e+13 GeV (metric)"),
 ("from about $9\\times10^{13}$ GeV for a constant rate to $3.1\\times10^{13}$ and $1.1\\times10^{14}$ GeV", "v6_quartic_decay.py", "plasma makes half at 3.13e+13 GeV (metric)"),
 ("from about $9\\times10^{13}$ GeV for a constant rate to $3.1\\times10^{13}$ and $1.1\\times10^{14}$ GeV", "v6_quartic_decay.py", "plasma makes half at 1.13e+14 GeV (metric)"),
 ("from $1.05\\times10^{14}$ to $5.2\\times10^{13}$ GeV", "fig_inflaton.py", "plasma makes half at 5.24e+13 GeV"),
 ("$C_M=I/\\pi^{3/2}=0.0115$", "v6_state_at_bang.py", "C_M = I/pi^1.5 = 0.0115"),
 ("ours cut off near $x=31$", "v6_hb_check.py", "cut off at x = 30.9 (Beltran-Palau et al. 0.06305)"),
 ("still carries 2 per cent", "v6_state_at_bang.py", "I/I_BNNP(0.06305) = 1.018"),
 ("$C_M=0.0023$", "v6_hb_check.py", "C_M = I_*/pi^1.5 = 0.0023"),
 ("$\\vartheta(1)=0.47$", "v6_hb_check.py", "Theta(1) = 0.471"),
 ("$c_H<4.1\\times10^{-5}$ at", "v6_bath.py", "M = 449 PeV: c < 4.1e-05"),
 ("$c_H<2.6\\times10^{-5}$ at $10^{14}$ GeV", "v6_bath.py", "M = 220 PeV: c < 2.6e-05"),
 ("under a tenth of a per cent of the dark matter", "v6_bath.py", "(0.062% of the dark matter)"),
 ("makes 40 times the dark matter", "v6_bath.py", "inflaton decays alone make 40 times"),
 ("a branching ratio above $3.5\\times10^{-20}$", "v6_bath.py", "branching to N1 pairs below 3.5e-20"),
 ("$y<9.2\\times10^{-12}\\,(T_{\\rm rh}/1.87\\times10^{13}\\,{\\rm GeV})^{1/2}", "v6_bath.py", "y_max = 9.2e-12 (T_rh/1.87e13 GeV)^(1/2) (449 PeV/M)^(1/2), whatever m_phi is"),
 ("is $0.05\\,M/M_P$ at the reference point", "v6_bath.py", "so y below 0.05 M/M_P keeps it under a tenth"),
 ("$r=1.6\\times10^{-5}$ at $H_{\\rm inf}=10^{12}$ GeV, where", "v6_bath.py", "r = 1.6e-05"),
 ("would carry $10^{78}$ times the plasma's energy", "v6_bath.py", "carries 1e+78 times"),
 ("from $4.2\\times10^5$ GeV up to about $H_{\\rm inf}$", "v6_bath.py", "the plasma sets M = 4.2e+05 GeV"),
 ("tail carries 9.1 per cent", "v6_reheating_production.py", "tail carries 9.1% of I"),
 ("The thermal mass exceeds the Hubble rate below $1.0\\times10^{16}$ GeV", "v6_bath.py", "that thermal mass exceeds H below T = 1.0e+16 GeV"),
 ("$(y\\phi)^4/16\\pi^2\\approx6.9\\times10^{61}$ GeV$^4$, would then exceed the plateau's $1.8\\times10^{61}$", "v6_bath.py", "one-loop energy 6.9e+61 GeV^4 against the plateau's 1.8e+61"),
 ("reheats perturbatively to at most about $3\\times10^{13}$ GeV", "v6_bath.py", "reheats to at most 3e+13 GeV"),
 ("the plasma shifts it by 0.4 per cent of $M_P$ in total", "v6_trilinear_bath.py", "displacement -3.84e-03 M_P (0.4 per cent)"),
 ("the plasma shifts it by 0.4 per cent of $M_P$ in total", "v6_trilinear_bath.py", "the plasma moves the inflaton by less than 1 per cent of M_P                   yes"),
 (r"leaves the horizon at $\phi=0.86\,M_P$", "v6_alpha_inflaton.py", "today's Hubble scale leaves at phi = 0.86 M_P"),
 ("suppressed on the plateau if reheating exceeds about $3\\times10^{13}$ GeV", "v6_bath.py", "reheats to at most 3e+13 GeV"),
 ("moves by 0.5 per cent for $H_{\\rm inf}=10^{11}$ GeV", "v6_start_check.py", "with no condition on m' the mass moves by -0.50 per cent"),
 ("changes $I$ by under $10^{-6}$", "v6_start_check.py", "the solver's start changes I by 2.2e-07"),
 ("plasma yield is 0.84 of the sudden formula", "v6_mass_smooth.py", "yield / sudden formula at T_max = 0.84"),
 ("gives a mass up to 44 per cent higher", "fig_inflaton.py", "up to 44.2 per cent higher"),
 ("between $1.3\\times10^{-2}$ and $4.9\\times10^{-2}$", "v6_state_at_bang.py", "Int x^3 p = 1.3e-02"),
 ("between $1.3\\times10^{-2}$ and $4.9\\times10^{-2}$", "v6_state_at_bang.py", "Int x^3 p = 4.9e-02"),
 ("$1.0\\times10^{-4}$ to $4.1\\times10^{-4}$ to change production", "v6_state_at_bang.py", "Int x^3 p = 1.0e-04 -> N = 13.9"),
 ("$1.0\\times10^{-4}$ to $4.1\\times10^{-4}$ to change production", "v6_state_at_bang.py", "6.6e-03, 4.1e-04"),
 ("against $1.6\\times10^{-3}$ to $6.6\\times10^{-3}$", "v6_state_at_bang.py", "1.4e-02, 1.6e-03"),
 ("it is 0.1 to 1.8 per cent over", "v6_state_at_bang.py", "the largest allowed change is 0.1%"),
 ("it is 0.1 to 1.8 per cent over", "v6_state_at_bang.py", "the largest allowed change is 1.8%"),
 ("11.5 to 12.9 for the state", "v6_state_at_bang.py", "N = 11.5 (bound)"),
 ("11.5 to 12.9 for the state", "v6_state_at_bang.py", "N = 12.9 (bound)"),
 ("12.0 to 13.4 for random phases", "v6_state_at_bang.py", "12.0 (up)"),
 ("12.0 to 13.4 for random phases", "v6_state_at_bang.py", "13.4 (up)"),
 ("12.7 to 14.1 for tuned phases", "v6_state_at_bang.py", "11.5 (bound), 12.0 (up), 12.7 (down)"),
 ("12.7 to 14.1 for tuned phases", "v6_state_at_bang.py", "12.9 (bound), 13.4 (up), 14.1 (down)"),
 ("| $10^{10}$ | $1.87\\times10^{13}$ | 12.7 | 13.2 | 13.9 | 13.8 |", "v6_state_at_bang.py", "N = 12.7 (bound), 13.2 (up), 13.8 (down)"),
 ("| $10^{10}$ | $1.87\\times10^{13}$ | 12.7 | 13.2 | 13.9 | 13.8 |", "v6_state_at_bang.py", "1.0e-04 -> N = 13.9"),
 ("| $10^{12}$ | $1.87\\times10^{13}$ | 11.9 | 12.3 | 13.1 | 13.1 |", "v6_state_at_bang.py", "N = 11.9 (bound), 12.3 (up), 13.1 (down)"),
 ("| $10^{12}$ | $1.87\\times10^{13}$ | 11.9 | 12.3 | 13.1 | 13.1 |", "v6_state_at_bang.py", "1.1e-04 -> N = 13.1"),
 ("| $10^{13}$ | $1.87\\times10^{13}$ | 11.5 | 12.0 | 12.7 | 12.7 |", "v6_state_at_bang.py", "N = 11.5 (bound), 12.0 (up), 12.7 (down)"),
 ("| $10^{13}$ | $1.87\\times10^{13}$ | 11.5 | 12.0 | 12.7 | 12.7 |", "v6_state_at_bang.py", "1.1e-04 -> N = 12.7"),
 ("| $10^{10}$ | $5\\times10^{13}$ | 12.9 | 13.4 | 14.1 | 14.1 |", "v6_state_at_bang.py", "N = 12.9 (bound), 13.4 (up), 14.1 (down)"),
 ("| $10^{10}$ | $5\\times10^{13}$ | 12.9 | 13.4 | 14.1 | 14.1 |", "v6_state_at_bang.py", "3.1e-04 -> N = 14.1"),
 ("| $10^{12}$ | $5\\times10^{13}$ | 12.1 | 12.6 | 13.3 | 13.3 |", "v6_state_at_bang.py", "N = 12.1 (bound), 12.6 (up), 13.3 (down)"),
 ("| $10^{12}$ | $5\\times10^{13}$ | 12.1 | 12.6 | 13.3 | 13.3 |", "v6_state_at_bang.py", "3.4e-04 -> N = 13.3"),
 ("| $10^{13}$ | $5\\times10^{13}$ | 11.7 | 12.2 | 13.0 | 12.9 |", "v6_state_at_bang.py", "N = 11.7 (bound), 12.2 (up), 12.9 (down)"),
 ("| $10^{13}$ | $5\\times10^{13}$ | 11.7 | 12.2 | 13.0 | 12.9 |", "v6_state_at_bang.py", "3.4e-04 -> N = 13.0"),
]

def check_table(note, fails):
    """Table 1 cell by cell against mass_table.py; masses are integers in PeV and the plasma share has
    one decimal."""
    for t, writer in TSVS.items():
        if not current(os.path.join(CALC, t), sources(writer), REC):
            fails.append("%s is older than %s or a module it imports" % (t, writer))
    import mass_table
    rows = re.findall(r"^\| \$([0-9.]+)\\times10\^\{(\d+)\}\$ \| (\d+) \| (\d+) \| ([0-9.]+)% \| (\d+) \| (\d+) \|$", note, re.M)
    if len(rows) < 10: fails.append("Table 1: only %d rows parsed" % len(rows)); return
    want = mass_table.rows()
    for a, e, *cells in rows:
        T = float(a)*10**int(e)
        r = [x for x in want if abs(x[0]/T - 1) < 1e-3]
        if not r: fails.append("Table 1: no solve at T_rh = %g" % T); continue
        for w, g, lab in zip(r[0][1:], cells, ("expansion", "with plasma", "share", "quartic", "sudden")):
            ok = ("%.1f" % w == g) if lab == "share" else ("%.0f" % w == g)
            if not ok: fails.append("Table 1 at %.3g GeV, %s: paper %s, script %.2f" % (T, lab, g, w))

def check(note, claims):
    fails, skipped = [], []
    body = flat(note)
    for ptxt, src, otxt in claims:
        if flat(ptxt) not in body: fails.append("note lacks %r" % ptxt)
        if src.startswith("log:"):
            path = os.path.join(CLOCK, src[4:])
            if not os.path.exists(path): skipped.append("%s (%s, the unpublished clock test)" % (ptxt, src[4:])); continue
            out = io.open(path, encoding="utf-8").read()
        else:
            here = os.path.exists(os.path.join(HERE, src.split()[0]))
            out = run(src, where=HERE if here else None)
        if otxt not in out: fails.append("%s output lacks %r (note: %r)" % (src, otxt, ptxt))
    check_table(note, fails)
    return fails, skipped

if __name__ == "__main__":
    note = io.open(NOTE, encoding="utf-8").read()
    if "--selftest" in sys.argv:
        bad = note.replace(r"$\alpha=0.0045$", r"$\alpha=0.0046$").replace("561 PeV or heavier", "600 PeV or heavier")
        bad = re.sub(r"^(\| \$3\\times10\^\{12\}\$ \| )(\d+)", lambda m: m.group(1) + str(int(m.group(2)) + 1), bad, count=1, flags=re.M)
        f, _ = check(bad, [c for c in CLAIMS if c[1] == A])
        ok = any("0.0045" in x for x in f) and any("561 PeV or heavier" in x for x in f) and any("Table 1" in x for x in f)
        print("selftest: a changed number in the text and a changed table cell are caught:", "yes" if ok else "NO")
        f, _ = check(note, [(r"Section 4 shows that couplings through the inflaton itself fail", A, "no such line in the output")])
        ok2 = len(f) == 1
        print("selftest: an output line the script does not print is caught:", "yes" if ok2 else "NO")
        sys.exit(0 if ok and ok2 else 1)
    fails, skipped = check(note, CLAIMS)
    for s in skipped: print("   skipped, not in this tree:", s)
    for f in fails: print("  ", f)
    print("inflaton paper: %d claims and Table 1 checked, %d skipped, %d problem(s)" % (len(CLAIMS) - len(skipped), len(skipped), len(fails)))
    sys.exit(1 if fails else 0)

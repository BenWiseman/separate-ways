#!/usr/bin/env python3
"""Table 1 of "Upside Down Inflation", printed as Markdown from the mass tables the solvers write:
checks/calc/v6_mass_smooth.tsv (perturbative reheating, Gamma = H(T_rh)) and
checks/calc/v6_mass_with_bath.tsv (sudden end of reheating) and v6_mass_smooth_w13.tsv (quartic
minimum, radiation-like reheating), so no cell is copied by hand. Masses are rounded to 1 PeV, the
precision the reheating modelling supports.

    python3 pub/inflaton/mass_table.py            # print the rows
"""
import os

HERE = os.path.dirname(os.path.abspath(__file__))
CALC = os.path.join(HERE, "..", "..", "checks", "calc")
ROWS = [3e12, 8e12, 1.2e13, 1.87e13, 2.5e13, 3.5e13, 5e13, 7e13, 1e14, 1.4e14, 2e14]

def tsv(name):
    out = {}
    with open(os.path.join(CALC, name)) as f:
        next(f)
        for line in f:
            p = line.rstrip("\n").split("\t")
            out.setdefault(p[0], {})[float(p[1])] = (float(p[3])/1e6, 100*float(p[4]))
    return out

def near(d, T):
    k = min(d, key=lambda t: abs(t/T - 1))
    if abs(k/T - 1) > 1e-3: raise SystemExit("no solve at T_rh = %g" % T)
    return d[k]

def latex_T(T):
    m, e = ("%.2e" % T).split("e")
    m = m.rstrip("0").rstrip(".")
    return "$%s\\times10^{%d}$" % (m, int(e))

def rows():
    p, s, q = tsv("v6_mass_smooth.tsv"), tsv("v6_mass_with_bath.tsv"), tsv("v6_mass_smooth_w13.tsv")
    out = []
    for T in ROWS:
        e = near(p["expansion only"], T)[0]
        w, share = near(p["with plasma"], T)
        su = near(s["central: alpha 1.2e-2, no boost"], T)[0]
        qu = near(q["with plasma"], T)[0]
        out.append((T, e, w, share, qu, su))
    return out

if __name__ == "__main__":
    print("| $T_{\\rm rh}$ (GeV) | expansion (PeV) | with plasma (PeV) | plasma share | quartic (PeV) | sudden end (PeV) |")
    print("|---|---|---|---|---|---|")
    for T, e, w, share, qu, su in rows():
        print("| %s | %.0f | %.0f | %.1f%% | %.0f | %.0f |" % (latex_T(T), e, w, share, qu, su))

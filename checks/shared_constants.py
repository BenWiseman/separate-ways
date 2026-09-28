#!/usr/bin/env python3
"""A number one script computes and others type must not drift from it.

    python3 checks/shared_constants.py
    python3 checks/shared_constants.py --selftest

WHY. J2end = 47.561945 is the non-degenerate Jacobi field at the far end of the contact geodesic.
contact_vanvleck.R computes it. Eight other files carry it as a literal, and it sits under two
things that matter: the calibration target Delta^{1/2} -> 3.9004 M s^{-1/2} that fork 9 has to
hit, and kappa = 0.0039, which sets the 2.1 micron shell. Eight copies of one computed number is
eight chances to go stale, and nothing was comparing them.

The value itself is not in question: vanvleck_by_a_second_route.R reproduces it to ten figures
from a closed-form geodesic and the analytic tidal field 3 M L^2 / r^5, sharing no machinery with
the numerically differentiated Riemann tensor that produced it. What this guards is the copies.
"""
import glob, io, re, subprocess, sys, os

SOURCE = "checks/calc/contact_vanvleck.R"
PATTERN = re.compile(r"defocusing:\s*J2\(pi\)\s*=\s*([0-9.]+)")
LITERAL = re.compile(r"\b47\.5\d*\b")


def computed():
    out = subprocess.run(["Rscript", SOURCE], capture_output=True, text=True,
                         cwd=os.getcwd()).stdout
    m = PATTERN.search(out)
    if not m:
        print("   %s no longer prints the value this gate reads  <-- ISSUE" % SOURCE)
        return None
    return float(m.group(1))


if "--selftest" in sys.argv:
    assert PATTERN.search("   defocusing: J2(pi) = 47.561945\n").group(1) == "47.561945", \
        "plant: the reader does not find the value in the line it is written on"
    assert not PATTERN.search("   defocusing: J2 = 47.561945"), \
        "plant: the reader matches a line it should not"
    assert LITERAL.findall("lam <- 3; J2end <- 47.561945  # x") == ["47.561945"], \
        "plant: the literal scan misses a copy"
    assert LITERAL.findall("nothing here") == [], "plant: the literal scan invents one"
    print("   plant: the reader and the literal scan both behave: yes")
    sys.exit(0)

val = computed()
if val is None:
    sys.exit(1)
bad, n = [], 0
for f in sorted(glob.glob("checks/calc/*.R") + glob.glob("checks/fig_*.R")):
    if f.endswith("contact_vanvleck.R"):
        continue
    for lit in LITERAL.findall(io.open(f, encoding="utf-8").read()):
        n += 1
        # the copy must be the computed value rounded to however many digits it carries
        dp = len(lit.split(".")[1])
        if abs(float(lit) - round(val, dp)) > 0.5 * 10 ** (-dp):
            bad.append((f.split("/")[-1], lit))
for f, lit in bad:
    print("   %s carries %s; %s computes %.6f  <-- ISSUE" % (f, lit, SOURCE.split("/")[-1], val))
if not bad:
    print("   J2end = %.6f, computed once and copied %d times, every copy agreeing" % (val, n))
sys.exit(1 if bad else 0)

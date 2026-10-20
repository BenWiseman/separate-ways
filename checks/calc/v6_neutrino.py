#!/usr/bin/env python3
"""Separate Ways v6, section 6: the neutrino masses the stabilising Z2 implies, and their tests.

With N1 removed from the seesaw the light-neutrino mass matrix has rank two and the lightest
neutrino is massless (Boyle, Finn and Turok). Inputs: NuFIT 6.0 with SK atmospheric data
(arXiv:2410.05380, Table 1): dm21 = 7.49(19)e-5, dm31 = +2.513(+21,-19)e-3 (NO),
dm32 = -2.484(20)e-3 (IO), s12^2 = 0.308, s13^2 = 0.02215 (NO, with SK).
"""
import math, sys
fails = 0
def check(ok, what):
    global fails
    print("   %-66s %s" % (what, "yes" if ok else "NO  <-- FAILED")); fails += (not ok)
dm21, e21 = 7.49e-5, 0.19e-5
dm31, e31p, e31m = 2.513e-3, 0.021e-3, 0.019e-3
dm32, e32 = 2.484e-3, 0.020e-3
S12, S13 = 0.308, 0.02215     # NuFIT 6.0 with SK, NO (Table 1)

def no_sum(a=dm21, b=dm31): return math.sqrt(a) + math.sqrt(b)
def io_sum(a=dm21, b=dm32): return math.sqrt(b) + math.sqrt(b - a)

s_no = no_sum()*1e3; s_io = io_sum()*1e3
up = math.hypot(no_sum(dm21+e21)-no_sum(), no_sum(dm21, dm31+e31p)-no_sum())*1e3
dn = math.hypot(no_sum()-no_sum(dm21-e21), no_sum()-no_sum(dm21, dm31-e31m))*1e3
eio = math.hypot(io_sum(dm21, dm32+e32)-io_sum(), io_sum(dm21+e21)-io_sum())*1e3
print("inputs: dm21 = %.2e, dm31 = %.3e (NO), dm32 = -%.3e (IO) eV^2" % (dm21, dm31, dm32))
print("normal ordering, m1 = 0:   Sum = %.2f +%.2f -%.2f meV" % (s_no, up, dn))
print("inverted ordering, m3 = 0: Sum = %.2f +- %.2f meV" % (s_io, eio))
c13 = 1 - S13
a = math.sqrt(dm21)*S12*c13*1e3; b = math.sqrt(dm31)*S13*1e3
print("0nubb, NO with m1 = 0: m_bb between %.2f and %.2f meV (Majorana phases free)" % (abs(a-b), a+b))
# inverted ordering with m3 = 0: m1 = sqrt(|dm32| - dm21), m2 = sqrt(|dm32|); one phase matters
m1i, m2i = math.sqrt(dm32 - dm21)*1e3, math.sqrt(dm32)*1e3
lo_io = abs(m1i*(1 - S12) - m2i*S12)*c13; hi_io = (m1i*(1 - S12) + m2i*S12)*c13
print("0nubb, IO with m3 = 0: m_bb between %.1f and %.1f meV (same inputs; s13^2 moves it by under 0.1 meV)" % (lo_io, hi_io))
# NuFIT 6.0 without SK atmospheric data (IC19): dm31 = 2.534e-3 for NO
print("normal ordering without SK data: Sum = %.2f meV" % (no_sum(dm21, 2.534e-3)*1e3))
check(abs(s_no - 58.78) < 0.02, "NO sum reproduces 58.78 meV")
check(abs(s_io - 98.92) < 0.03, "IO sum reproduces 98.92 meV")
check(abs((a+b) - 3.72) < 0.05 and abs(abs(a-b) - 1.50) < 0.05, "m_bb band reproduces 1.50 to 3.72 meV")
# plant: an IO computed with the NO splitting must differ (the v5.6 error)
check(abs(io_sum(dm21, dm31)*1e3 - s_io) > 0.2, "plant: using the NO splitting for IO gives a different number")
sys.exit(1 if fails else 0)

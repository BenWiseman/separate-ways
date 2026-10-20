#!/usr/bin/env python3
"""Every number Separate Ways v6 quotes must be in the output of the script that makes it.

    python3 papers/1_separate_ways_v6/v6_numbers_check.py              # check (script outputs cached by mtime)
    python3 papers/1_separate_ways_v6/v6_numbers_check.py --selftest   # planted wrong numbers must fail

Each claim is (text in the paper, script, text in that script's output); the rows of Table 1 and
Table A.1 are claims like any other. Literature values (NuFIT, DESI, the KM3NeT event, the lifetime
bounds, the published rates of Appendix C) are checked against their sources elsewhere and are not
listed here. The second paper's check, pub/inflaton/inflaton_check.py, runs its scripts through run()
below.
"""
import glob, hashlib, io, os, re, subprocess, sys, tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
CALC = os.path.join(ROOT, "checks", "calc")
PAPER = os.path.join(HERE, "SEPARATE_WAYS_v6.md")
PY = os.path.join(ROOT, ".venv", "bin", "python3")
if not os.path.exists(PY): PY = sys.executable        # the release repository has no .venv
CACHE = os.path.join(os.environ.get("TMPDIR", "/tmp"), "v6_numbers_cache")

CLAIMS = [
 # abstract
 (r"$1.20\times10^{-2}\,T^8/M_P^4$ for a Majorana fermion", "v6_grav_rate.py", "alpha = 4 x + 22.5 x + 12 x = 1.2027e-02"),
 ("which by itself makes the coefficients 1.41\ntimes too large for fermion dark matter and 1.38 times", "v6_grav_rate.py", "A12): 1.6973e-02, 1.41 times ours"),
 ("which by itself makes the coefficients 1.41\ntimes too large for fermion dark matter and 1.38 times", "v6_grav_rate.py", "with 45 Dirac 2.9883e-03, 1.381 times too high"),
 ("rebuilt from their equations, is 4.27 times ours", "v6_grav_rate.py", "(4.27 x ours)"),
 (r"ever reached $8.0\times10^{13}$ GeV makes all the dark matter", "v6_bath.py", "= 8.008e+13 GeV"),
 ("491.6 PeV bounds the mass from\nabove", "v6_abundance.py", "M_* = 491.6 PeV"),
 ("its neutrino line lies at or below 245.8 PeV; with", "v6_abundance.py", "line at M_*/2 = 245.8 PeV"),
 ("the neutrino masses sum to 58.8 meV in normal ordering", "v6_neutrino.py", "Sum = 58.78 +0.24 -0.22 meV"),
 # sections 2 and 3
 ("$I_*=0.0127597$", "v6_no_floor.py", "I_* = 0.0127596674"),
 ("$I_{\\rm in}=1/(4\\pi^{3/2})=0.04490$", "v6_reheating_production.py", "0.044897"),
 ("and $I=0.0642$", "v6_state_at_bang.py", "I = 0.0642, C_M"),
 (r"M_*=491.6\ {\rm PeV}", "v6_abundance.py", "M_* = 491.6 PeV"),
 ("corresponds to $\\Omega_ch^2=0.1198$", "v6_abundance.py", "Omega_c h^2 = 0.11981"),
 ("gives $491.9\\pm2.0$ PeV", "v6_abundance.py", "= 491.9 PeV"),
 ("gives $491.9\\pm2.0$ PeV", "v6_abundance.py", "+- 2.0 PeV (mass)"),
 ("which on its own gives 484.8 PeV", "v6_abundance.py", "with s0 = 2.3e-38: 484.8 PeV"),
 ("$F=2$ gives 648.7 PeV", "v6_abundance.py", "F = 2: 648.7 PeV"),
 ("the ceiling rises by 1.6 PeV", "v6_abundance.py", "g* = 110.25: +1.6 PeV"),
 ("330.4 PeV at $x_c=0.564$", "v6_no_floor.py", "0.564            1.00        2.700    330.4"),
 ("gives the ceiling, 491.6 PeV", "v6_abundance.py", "M_* = 491.6 PeV"),
 ("which gives 330.4 PeV (Table A.1)", "v6_no_floor.py", "0.564            1.00        2.700    330.4"),
 ("89.6 PeV at $x_c=1.693$", "v6_no_floor.py", "1.693            3.00       70.467     89.6"),
 # section 4
 ("\\frac{N_{\\rm s}}{648000}+\\frac{49\\,N_{\\rm f}}{13824000}+\\frac{N_{\\rm v}}{54000}\\Bigr)=1.20\\times10^{-2}", "v6_grav_rate.py", "Dirac fermion          r = 1/648000, 49/13824000, 1/54000; R = 1.2027e-02"),
 ("$11351\\pi^3/20736000=1.697\\times10^{-2}$ per event, 1.41 times $\\alpha$, and Clery et al.'s Eq. (A12) is exactly twice that", "v6_grav_rate.py", "11351 pi^3/20736000 = 1.697301e-02"),
 ("Clery et al.'s Eq. (A12) is exactly twice that", "v6_grav_rate.py", "3.394601e-02, 2.000 x our 45-Dirac event sum"),
 ("twice our 45-Dirac sum for scalar dark matter; that sum, $3997\\pi^3/41472000$, is 1.38 times the corrected value", "v6_grav_rate.py", "3997 pi^3/20736000 = 5.976663e-03, 2.000 x our 45-Dirac scalar event sum"),
 ("by $1.41^{1/3}=1.12$ for fermion dark matter", "v6_grav_rate.py", "1.41^(1/3) = 1.122"),
 ("is 4.27 times $\\alpha$", "v6_grav_rate.py", "(4.27 x ours)"),
 ("gives a factor 1.88", "v6_grav_rate.py", "vertex 1.882"),
 ("the remaining 1.135 is 1.018 from the missing", "v6_grav_rate.py", "scalar factor and statistics 1.135"),
 ("coefficients to $1.0\\times10^{-4}$; it keeps", "v6_grav_rate.py", "so their bound f2 < 7e-5 rises to 9.99e-05"),
 ("lower $C_{\\rm UV}$ by 2.56", "v6_grav_rate.py", "their coefficients over ours: Einstein limit 4.272, above T* 2.558"),
 ("falls under one per cent for $f_2<2.15\\times10^{-5}$", "v6_grav_rate.py", "under one per cent of the dark matter for f2 < 2.15e-05"),
 ("$\\Delta m^2_{21}=7.49\\times10^{-5}$ eV$^2$ and $\\Delta m^2_{31}=2.513\\times10^{-3}$ eV$^2$", "v6_neutrino.py", "inputs: dm21 = 7.49e-05, dm31 = 2.513e-03 (NO), dm32 = -2.484e-03 (IO)"),
 ("| one Majorana $N_1$ per Dirac-pair event | $1.37\\times10^{-2}$ | 1.135 |", "v6_figures.py", "BTV steps in units of alpha 4.273, 2.270, 1.135, 1.115, 1.099, 1.000"),
 ("| both in Maxwell-Boltzmann statistics | $1.1805$ against $1.1560\\times10^{-2}$ | 1.021 |", "v6_grav_rate.py", "Maxwell-Boltzmann throughout: ours 1.1560e-02, corrected BTV 1.1805e-02 (ratio 1.021)"),
 ("| real scalar, $\\xi=0$ | $\\pi^3/648000$ | $49\\pi^3/82944000$ | $\\pi^3/324000$ | $2.16\\times10^{-3}$ | 1.38 |", "v6_grav_rate.py", "real scalar, xi = 0    r = 1/648000, 49/82944000, 1/324000; R = 2.1641e-03; particles per event 2; with 45 Dirac 2.9883e-03, 1.381"),
 ("| real scalar, $\\xi=1/6$ | $\\pi^3/3888000$ | $49\\pi^3/82944000$ | $\\pi^3/324000$ | $2.00\\times10^{-3}$ | 1.41 |", "v6_grav_rate.py", "real scalar, xi = 1/6  r = 1/3888000, 49/82944000, 1/324000; R = 2.0046e-03; particles per event 2; with 45 Dirac 2.8288e-03, 1.411"),
 ("| Dirac fermion | $\\pi^3/648000$ | $49\\pi^3/13824000$ | $\\pi^3/54000$ | $1.20\\times10^{-2}$ | 1.41 |", "v6_grav_rate.py", "Dirac fermion          r = 1/648000, 49/13824000, 1/54000; R = 1.2027e-02; particles per event 2; with 45 Dirac 1.6973e-02, 1.411"),
 ("| massless vector | $\\pi^3/324000$ | $49\\pi^3/6912000$ | $\\pi^3/27000$ | $2.41\\times10^{-2}$ | 1.41 |", "v6_grav_rate.py", "massless vector        r = 1/324000, 49/6912000, 1/27000; R = 2.4055e-02; particles per event 2; with 45 Dirac 3.3946e-02, 1.411"),
 ("enters only the channel with plasma scalars", "v6_grav_rate.py", "only the scalar channel depends on the dark matter's xi (fermion and gauge stress tensors are traceless) yes"),
 ("=480.6", "v6_bath.py", "C = 480.6"),
 ("=8.9\\times10^{-19}", "v6_bath.py", "Y_obs(M_c) = 8.882e-19"),
 (r"T_\star=8.0\times10^{13}", "v6_bath.py", "= 8.008e+13 GeV"),
 ("a temperature 4.3 times the $1.87\\times10^{13}$ GeV at which $H=M_*$", "v6_bath.py", "T*/T(H=M_c) = 4.28"),
 ("$1.95\\times10^3$ times the observed", "v6_bath.py", "1.95e+03 times the dark matter"),
 ("$2.8\\times10^{13}$ times at $T=M_P$", "v6_bath.py", "2.81e+13 times the dark matter"),
 ("$\\Gamma=1.9\\times10^{-2}H$", "v6_bath.py", "T = 2.44e+18: Gamma/H = 1.9e-02"),
 ("252 GeV if the plasma", "v6_bath.py", "M = 252 GeV"),
 ("$7.42\\times10^{-2}\\,T^8/M_P^4$, 6.2 times", "v6_grav_rate.py", "alpha_EC = 7.4185e-02 = 6.17 alpha"),
 ("lowers $T_\\star$ to $4.2\\times10^{13}$ GeV", "v6_grav_rate.py", "from 8.008e13 to 4.15e+13 GeV"),
 ("$7.4\\times10^{-6}\\,c_H^2", "v6_bath.py", "Y_portal = 7.4e-06"),
 ("$c_H>3.5\\times10^{-7}$", "v6_bath.py", "c = 3.5e-07 makes all the dark matter"),
 ("about $2\\times10^{-10}$, and the operator never matters", "v6_bath.py", "c ~ M/M_P = 2.0e-10"),
 ("$\\alpha$ rises by 1.8 per cent", "v6_grav_rate.py", "alpha rises by 1.8 per cent"),
 # section 5
 ("weighs at most 491.6 PeV", "v6_abundance.py", "M_* = 491.6 PeV"),
 ("58.78^{+0.24}_{-0.22}", "v6_neutrino.py", "Sum = 58.78 +0.24 -0.22 meV"),
 ("98.92\\pm0.40", "v6_neutrino.py", "Sum = 98.92 +- 0.40 meV"),
 ("neutrino under $10^{-55}$ eV", "v6_line_rates.py", "< 2e-56 eV"),
 ("needs $y<6\\times10^{-31}$", "v6_line_rates.py", "needs y < 5.8e-31"),
 ("0.14 events in IceCube's 12.6-year", "v6_line_rates.py", "M = 492 PeV, tau = 1.00e+29 s: IceCube 12.6 yr 0.141 events"),
 ("and 0.64 in ten years", "v6_line_rates.py", "M = 492 PeV, tau = 1.00e+29 s: IceCube 12.6 yr 0.141 events, Gen2 optical 10 yr 0.64"),
 ("a cored halo lowers them by 15 per cent", "v6_line_rates.py", "0.85 of the NFW one"),
 ("Gen2 would see one\nevent at $6\\times10^{28}$ s", "v6_line_rates.py", "for tau = 6.4e+28 s"),
 ("at most 245.8 PeV,\nso under the conditions of Section 4", "v6_abundance.py", "line at M_*/2 = 245.8 PeV"),
 # section 6
 ("a factor 1.41 below what the 45 count gives", "v6_grav_rate.py", "A12): 1.6973e-02, 1.41 times ours"),
 ("is 4.27 times ours (Fig. 1b)", "v6_grav_rate.py", "(4.27 x ours)"),
 ("once its plasma has reached $8.0\\times10^{13}$ GeV", "v6_bath.py", "= 8.008e+13 GeV"),
 ("their mass is a ceiling, 491.6 PeV, with no floor", "v6_abundance.py", "M_* = 491.6 PeV"),
 ("well below $1.0\\times10^{-4}$", "v6_grav_rate.py", "rises to 9.99e-05"),
 ("its neutrino line lies at or below 245.8 PeV", "v6_abundance.py", "line at M_*/2 = 245.8 PeV"),
 ("sum to 58.8 meV in normal ordering. With inflation", "v6_neutrino.py", "Sum = 58.78 +0.24 -0.22 meV"),
 # appendix B and section 4
 ("so their thermal abundances are 1.41 times too high for Dirac fermion dark matter and 1.38 times", "v6_grav_rate.py", "published thermal abundances: high by 1.411 (fermion) and 1.381 (scalar)"),
 ("so in their rate those two channels stand at half of ours and the fermion channel at twice", "v6_grav_rate.py", "Garny et al.'s channels in their rate, relative to ours: scalar 0.50, fermion 2.00, vector 0.50"),
 ("The vector row is twice the Dirac row, as the widths", "v6_grav_rate.py", "massless vector        r = 1/324000, 49/6912000, 1/27000; R = 2.4055e-02"),
 ("The vector row is twice the Dirac row, as the widths", "v6_grav_rate.py", "with 1/2 for identical pairs the channels stand 1 : 6.000 : 12.000"),
 # appendix C
 ("In their rate the three channels therefore stand at 1/2, 2 and 1/2 of ours", "v6_grav_rate.py", "Garny et al.'s channels in their rate, relative to ours: scalar 0.50, fermion 2.00, vector 0.50"),
 ("a quarter of ours, ours, and a quarter of ours", "v6_grav_rate.py", "Garny et al.'s channels over ours (scalar, fermion, vector): 0.2500, 1.0000, 0.2500 at every angle"),
 ("is $1/16$ of our scalar channel", "v6_grav_rate.py", "Bernal et al.'s channels over ours (scalar, fermion, vector): 0.0625, 0.0625, 0.0625 at every angle"),
 ("are ours with the factor $1/2$ for an identical pair; the fermion channel", "v6_grav_rate.py", "Clery et al.'s channels over ours (scalar, fermion, vector): 0.5000, 1.0000, 0.5000 at every angle"),
 ("$11351\\pi^3/331776000=1.06\\times10^{-3}$, exactly $1/16$ of our 45-Dirac sum", "v6_grav_rate.py", "Bernal et al.'s Eq. (B17) is 11351 pi^3/331776000 = 1.0608e-03, 0.0882 x ours, 0.0625 of our 45-Dirac sum"),
 ("per event it is $1.697\\times10^{-2}$, 1.41 times $\\alpha$", "v6_grav_rate.py", "A12): 1.6973e-02, 1.41 times ours"),
 ("$s^2/480$, $s^2/160$ and $s^2/40$", "v6_grav_rate.py", "in Tang and Wu's normalisation our channels are A = s^2/480, s^2/160, s^2/40"),
 ("$\\pi/30:\\pi/5:2\\pi/5=1:6:12$", "v6_grav_rate.py", "with 1/2 for identical pairs the channels stand 1 : 6.000 : 12.000"),
 ("carries half its sum (`v6_grav_rate.py`)", "v6_grav_rate.py", "a Weyl pair (one helicity combination) carries exactly half the Dirac sum yes"),
 ("| 1/16 each | 1/16, 1/8, 1/16 | 45 Dirac | $1.06\\times10^{-3}$ | 0.088 |", "v6_grav_rate.py", "1.0608e-03, 0.0882 x ours"),
 ("| 1, 2, 1 | 45 Dirac | $3.39\\times10^{-2}$ | 1.41 |", "v6_grav_rate.py", "11351 pi^3/10368000 = 3.394601e-02, 2.000 x our 45-Dirac event sum"),
 ("| 1, 1, 29/12 | Table B.1 | 47 Weyl | $5.14\\times10^{-2}$ | 2.14 |", "v6_grav_rate.py", "per event (one pair per collision): 2.136 x ours"),
 ("| this paper | 1, 1, 1 | 1, 1, 1 | 22.5 Dirac | $1.20\\times10^{-2}$ | 1 |", "v6_grav_rate.py", "alpha = 4 x + 22.5 x + 12 x = 1.2027e-02"),
 # appendix B
 ("=1.2027\\times10^{-2}", "v6_grav_rate.py", "alpha = 4 x + 22.5 x + 12 x = 1.2027e-02"),
 ("$11351\\pi^3/20736000=1.697\\times10^{-2}$ per event, exactly half of Eq. (A12)", "v6_grav_rate.py", "3.394601e-02, 2.000 x our 45-Dirac event sum"),
 ("$C_{\\rm IR}\\approx2.1\\times10^{-3}$", "v6_grav_rate.py", "C_IR = 2.142e-03"),
 ("=1.78\\times10^{-3}$", "v6_grav_rate.py", "= 1.782e-03"),
 ("4.27 times ours", "v6_grav_rate.py", "(4.27 x ours)"),
 ("their rate is 1.115 times", "v6_grav_rate.py", "(1.115 x ours)"),
 ("the two agree to 2 per cent", "v6_grav_rate.py", "(ratio 1.021)"),
 ("their bath's 47 Weyl fermions", "v6_grav_rate.py", "the 2 per cent Boltzmann residual is BTV's 47 Weyl fermions against 45   yes"),
 ("| as printed | $5.14\\times10^{-2}$ | 4.27 |", "v6_grav_rate.py", "as printed: 5.1386e-02 (4.27 x ours)"),
 ("| transverse gauge-boson vertex | $2.73\\times10^{-2}$ | 2.27 |", "v6_grav_rate.py", "transverse gauge vertex: 2.7305e-02"),
 ("| $1.37\\times10^{-2}$ | 1.135 |", "v6_grav_rate.py", "Majorana pair counting: 1.3652e-02"),
 ("| $1.34\\times10^{-2}$ | 1.115 |", "v6_grav_rate.py", "identical real scalars: 1.3416e-02 (1.115 x ours)"),
 ("integrates over angle to $3\\pi$", "v6_grav_rate.py", "contact channel integrates to 3 pi                                       yes"),
 ("fifteen times the graviton's Dirac channel", "v6_grav_rate.py", "15.000 times the graviton's Dirac channel"),
 ("$1.9\\times10^{-2}H$ at $T=M_P$, grows as", "v6_grav_rate.py", "Gamma/H = 1.9e-02 at the reduced Planck temperature"),
 ("would pass it\nat $3.7\\,M_P$", "v6_grav_rate.py", "1 at 3.7 M_P"),
 ("at $5.3\\times10^{13}$ GeV for $f_2=1.0\\times10^{-4}$, where it is $1.9\\times10^{-16}H$", "v6_grav_rate.py", "T* = 5.254e+13 GeV, where Gamma/H = 1.9e-16"),
 ("$C_{\\rm UV}=C_{\\rm IR}/768$", "v6_grav_rate.py", "their C_UV is exactly C_IR/768"),
 ("They lower $C_{\\rm UV}$ by\n2.56 against 4.27 for $C_{\\rm IR}$", "v6_grav_rate.py", "their coefficients over ours: Einstein limit 4.272, above T* 2.558"),
 ("rises by $(4.27^{1/4}\\,2.56^{3/4})^{1/3}=1.43$, to $1.0\\times10^{-4}$", "v6_grav_rate.py", "by the factor (4.27^(1/4) 2.56^(3/4))^(1/3) = 1.43"),
 ("| as printed | 2.558 |", "v6_grav_rate.py", "above T* over ours, step by step:      2.558, 1.359, 0.680, 0.668, 0.658, 1.000"),
 ("| transverse gauge-boson vertex | 1.359 |", "v6_grav_rate.py", "above T* over ours, step by step:      2.558, 1.359"),
 ("| one Majorana $N_1$ per Dirac-pair event | 0.680 |", "v6_grav_rate.py", "2.558, 1.359, 0.680"),
 ("| $1/2$ for identical real scalars | 0.668 |", "v6_grav_rate.py", "2.558, 1.359, 0.680, 0.668"),
 ("the mean of $1/s$ over the rate's integrand is $0.051/T^2$", "v6_grav_rate.py", "<1/s> = 0.0514/T^2 for bosons"),
 ("$g_3^2=0.32$ at $10^{14}$ GeV", "v6_grav_rate.py", "at 1e14 GeV g3^2 = 0.319"),
 ("3.8 of it from errors", "v6_grav_rate.py", "the errors make 3.83"),
 ("| 45 Weyl fermions, not 47 | 0.658 |", "v6_grav_rate.py", "0.668, 0.658, 1.000"),
 ("1.018 from the missing $1/2$", "v6_grav_rate.py", "1.135 = 1.018 (1/2 for identical scalars)"),
 ("1.015 from their\nbath's 47 Weyl fermions", "v6_grav_rate.py", "x 1.015 (47 Weyl fermions)"),
 ("1.099 from their Maxwell-Boltzmann average", "v6_grav_rate.py", "x 1.099 (their statistics)"),
 ("give 1.015\nof that", "v6_grav_rate.py", "x 1.015 (47 Weyl fermions)"),
 ("$2.3\\times10^{-3}$ of the radiation's energy", "v6_no_floor.py", "holds 2.3e-03 of the radiation's energy"),
 ("$\\langle s\\rangle/m_D^2\\approx69$", "v6_grav_rate.py", "= 44/0.64 = 69"),
 ("still $2\\times10^{-14}$ per $N_1$", "v6_grav_rate.py", "so 2e-14 per N1 across the crossing"),
 ("add 1.5, 1.8 or 2.1 per cent", "v6_grav_rate.py", "raise the rate by 1.5 per cent in their statistics, 1.8 in ours"),
 ("add 1.5, 1.8 or 2.1 per cent", "v6_grav_rate.py", "(ratio 1.021)"),
 ("on $f_2$ by 0.12 per cent", "v6_grav_rate.py", "moves that bound by 0.12 per cent"),
 ("$2.28\\times10^{22}$ GeV cm$^{-2}$", "v6_line_rates.py", "that is 2.28e+22 GeV cm^-2"),
 ("which is $1393\\zeta(3)/(25216\\pi^3)$", "v6_grav_rate.py", "the same as their C_IR = 2.142e-03"),
 ("$2.5\\times10^{-16}$ interactions", "v6_grav_rate.py", "each N1 meets 2.5e-16 graviton-exchange interactions"),
 ("$x\\approx10^5$ for the ceiling's mass", "v6_grav_rate.py", "has x = 1.1e+05"),
 ("| 45 Weyl fermions, not 47 | $1.32\\times10^{-2}$ | 1.099 |", "v6_grav_rate.py", "with 45 Weyl fermions in place of their 47: 1.3217e-02 (1.099 x ours)"),
 ("with our gauge-boson channel it would be\n$1.1\\times10^{-3}$", "v6_grav_rate.py", "it would be 1.14e-03, against the 2.14e-03 they quote"),
 ("lower the\nrow by 0.6 per cent", "v6_grav_rate.py", "massive-vector row falls by 0.6 per cent"),
 ("lowers it threefold against the photons", "v6_no_floor.py", "against the photons' temperature by 3.01"),
 ("the whole rate moves at the order of 1 per cent", "v6_grav_rate.py", "the rate shifts at relative order 0.0102"),
 ("| massive vector, $T\\gg m$ | $\\pi^3/216000$ | $637\\pi^3/82944000$ | $13\\pi^3/324000$ | $2.62\\times10^{-2}$ | 1.41 |", "v6_grav_rate.py", "massive vector, T >> m r = 1/216000, 637/82944000, 13/324000; R = 2.6219e-02; particles per event 2; with 45 Dirac 3.6934e-02, 1.409 times too high"),
 ("Eq. (A8), $3997\\pi^3/20736000$", "v6_grav_rate.py", "Clery et al.'s published Eq. (A8) is 3997 pi^3/20736000"),
 ("gauge-boson channel, 29/12 of ours, gives a factor 1.88", "v6_grav_rate.py", "their (D4) is 29/12 of our gauge-boson channel"),
 ("equals the plasma temperature when\nthe mass is 0.26 keV", "v6_no_floor.py", "M = 0.26 keV"),
 ("moves the mass by 0.12 per cent", "v6_no_floor.py", "at w = 0.01 it is 0.12 per cent"),
 ("integrates to $2\\pi/15$, six times the conformally coupled value", "v6_grav_rate.py", "6.000 times xi = 1/6"),
 ("integrates to $29\\pi/15$, against our $4\\pi/5$", "v6_grav_rate.py", "integrates to 6.0737 = 29 pi/15"),
 ("physical polarisations give 4.17 times our channel, or 1.04", "v6_grav_rate.py", "physical polarisations 4.167, Feynman-gauge sum 6.833"),
 ("physical polarisations give 4.17 times our channel, or 1.04", "v6_grav_rate.py", "physical polarisations 1.042, Feynman-gauge sum 1.708"),
 ("Feynman-gauge sum 6.83 or 1.71", "v6_grav_rate.py", "Feynman-gauge sum 6.833"),
 ("Feynman-gauge sum 6.83 or 1.71", "v6_grav_rate.py", "Feynman-gauge sum 1.708"),
 ("give $43/(10\\pi^5)$, 1.17 times $\\alpha$", "v6_grav_rate.py", "43/(10 pi^5), 1.168 x ours"),
 ("1.12 times $\\alpha$ with quantum statistics", "v6_grav_rate.py", "1.117 x ours (they print no total)"),
 ("| 1/2, 2, 1/2 | 45 Dirac | none | 1.17 ($g=1$) |", "v6_grav_rate.py", "43/(10 pi^5), 1.168 x ours"),
 ("$f=0.116$ neutrinos at the full energy", "v6_line_rates.py", "(0.116 line neutrinos per decay)"),
 ("$D=7.4$ GeV cm$^{-3}$ kpc", "v6_line_rates.py", "sky-averaged column 7.38 GeV cm^-3 kpc"),
 ("$A=1.99\\times10^3$ m$^2$", "v6_line_rates.py", "1.99e+03 m^2 (IceCube EHE)"),
 ("$A=1.14\\times10^4$ m$^2$", "v6_line_rates.py", "1.14e+04 m^2 (Gen2 optical)"),
 ("which vanishes for one choice of them and grows without limit", "v6_grav_rate.py", "the contact term vanishes at alpha = 0, beta = +-1 and grows without limit yes"),
 ("gives 1.6 per cent in its\nchannel", "v6_grav_rate.py", "so the gluon's m^2 <1/s> = 0.016"),
 ("it gives $7.42\\times10^{-2}$, 6.17 times", "v6_grav_rate.py", "alpha_EC = 7.4185e-02 = 6.17 alpha"),
 ("3 per cent below what Shaposhnikov", "v6_grav_rate.py", "(1.031 of ours"),
 ("| 0.141 | 1.008 | 490.1 |", "v6_no_floor.py", "0.141            0.25        1.008    490.1     245.0"),
 ("| 0.282 | 1.119 | 470.0 |", "v6_no_floor.py", "0.282            0.50        1.119    470.0     235.0"),
 ("| 0.423 | 1.574 | 410.0 |", "v6_no_floor.py", "0.423            0.75        1.574    410.0     205.0"),
 ("| 0.564 | 2.700 | 330.4 |", "v6_no_floor.py", "0.564            1.00        2.700    330.4     165.2"),
 ("| 0.846 | 8.312 | 210.7 |", "v6_no_floor.py", "0.846            1.50        8.312    210.7     105.4"),
 ("| 1.128 | 20.256 | 147.6 |", "v6_no_floor.py", "1.128            2.00       20.256    147.6      73.8"),
 ("| 1.693 | 70.467 | 89.6 |", "v6_no_floor.py", "1.693            3.00       70.467     89.6      44.8"),
]

# Scripts that integrate modes take from minutes to hours, so their saved logs are read instead of
# rerunning them. A source may carry the flags it was run with ("v6_mass_smooth.py --w13" reads
# logs/v6_mass_smooth_w13.log). A log is refused if the script, or any v6 module it imports at any
# depth, changed after the log was made.
SLOW = {"v6_mass_smooth.py", "v6_reheating_smooth.py", "v6_hinf.py", "v6_mass_with_bath.py",
        "v6_state_at_bang.py", "v6_instant.py", "v6_reheating_production.py", "v6_start_check.py",
        "v6_quartic_decay.py"}
# each mass table and the script that writes it
TSVS = {"v6_mass_smooth.tsv": "v6_mass_smooth.py", "v6_mass_smooth_w13.tsv": "v6_mass_smooth.py",
        "v6_mass_smooth_ec.tsv": "v6_mass_smooth.py", "v6_mass_smooth_w13_ec.tsv": "v6_mass_smooth.py",
        "v6_mass_with_bath.tsv": "v6_mass_with_bath.py"}

def tree(script, seen=None):
    """the script and every v6 module it imports, at any depth"""
    seen = set() if seen is None else seen
    if script not in seen:
        seen.add(script)
        text = io.open(os.path.join(CALC, script), encoding="utf-8").read()
        for m in re.findall(r"(?m)^\s*(?:from|import)\s+(v6_[a-z0-9_]+)", text):
            tree(m + ".py", seen)
    return seen

# A git checkout stamps every file with the moment it was checked out, so in the release repository
# the ages of files say nothing about which was made first. The export writes logs/SOURCES.sha256, the
# hash of each source as exported, and a source that still has its recorded hash counts as unchanged
# since the logs were made, whatever its age. In the working tree there is no record and age decides.
HASHES = os.path.join(CALC, "logs", "SOURCES.sha256")
REC = {}
if os.path.exists(HASHES):
    REC = {n: h for h, n in (l.split() for l in io.open(HASHES, encoding="utf-8") if l.strip())}

def sha(path):
    return hashlib.sha256(io.open(path, "rb").read()).hexdigest()

def current(target, sources, rec):
    """a log or mass table is current if every source still has its recorded hash or is older than it"""
    if not os.path.exists(target): return False
    t = os.path.getmtime(target)
    return all(rec.get(os.path.basename(x)) == sha(x) or os.path.getmtime(x) <= t for x in sources)

def sources(script):
    return [os.path.join(CALC, x) for x in sorted(tree(script))]

def run(source, where=None):
    """the output of a script, from its saved log if it is slow; where names the directory of a script
    that sits beside its own paper, as the second paper's do"""
    script, *flags = source.split()
    if script in SLOW:
        name = script.replace(".py", "") + ("_w13" if "--w13" in flags else "") + ("_ec" if "--ec" in flags else "")
        log = os.path.join(CALC, "logs", name + ".log")
        if not current(log, sources(script), REC):
            raise SystemExit("%s: saved log missing, or the script or a module it imports changed after it was made;"
                             " rerun it into %s" % (source, log))
        out = io.open(log, encoding="utf-8").read()
        if "exit 0" not in out: raise SystemExit("%s: saved log does not end in exit 0" % source)
        return out
    # the figure script sits beside the paper; every other script is in calc
    if where: path = os.path.join(where, script)
    else: path = os.path.join(HERE, script) if script.startswith("v6_figures") else os.path.join(CALC, script)
    deps = [path] + sorted(glob.glob(os.path.join(CALC, "v6_*.py"))) + [os.path.join(CALC, t) for t in TSVS]
    key = hashlib.sha1("".join("%s%.0f" % (d, os.path.getmtime(d)) for d in sorted(deps)).encode()).hexdigest()[:12]
    os.makedirs(CACHE, exist_ok=True)
    cf = os.path.join(CACHE, "%s.%s.out" % (script, key))
    if os.path.exists(cf): return io.open(cf, encoding="utf-8").read()
    r = subprocess.run([PY, path], capture_output=True, text=True, cwd=ROOT, timeout=7200)
    if r.returncode != 0: raise SystemExit("%s exited %d:\n%s" % (script, r.returncode, r.stdout[-2000:] + r.stderr[-2000:]))
    io.open(cf, "w", encoding="utf-8").write(r.stdout); return r.stdout

def flat(t):
    """Whitespace runs to one space, so a reflowed line does not break a claim."""
    return " ".join(t.split())

def check(paper, claims):
    fails = []
    body = flat(paper)
    for ptxt, script, otxt in claims:
        if flat(ptxt) not in body: fails.append("paper lacks %r" % ptxt)
        if otxt is None: continue
        out = run(script)
        if otxt not in out: fails.append("%s output lacks %r (paper: %r)" % (script, otxt, ptxt))
    return fails

if __name__ == "__main__":
    paper = io.open(PAPER, encoding="utf-8").read()
    if "--selftest" in sys.argv:
        bad = paper.replace("T_\\star=8.0\\times10^{13}", "T_\\star=8.1\\times10^{13}")
        # plant a table error: one cell of Table A.1
        bad = bad.replace("| 0.564 | 2.700 | 330.4 |", "| 0.564 | 2.700 | 330.5 |")
        if bad.count("T_\\star=8.1") != 1 or "330.5 |" not in bad: sys.exit("selftest could not plant its errors")
        f = check(bad, [c for c in CLAIMS if c[1] in ("v6_bath.py", "v6_no_floor.py")])
        ok = any("8.0" in x for x in f) and any("330.4 |" in x for x in f)
        print("selftest: a changed number in the text and a changed table cell are both caught:", "yes" if ok else "NO")
        # the freshness rule: v6_hinf.py reaches v6_bath.py only through v6_mass_with_bath.py, and a log
        # older than its source, as every file can be after a checkout, is current only while the
        # source has the hash recorded for it
        deep = "v6_bath.py" in tree("v6_hinf.py")
        with tempfile.NamedTemporaryFile() as tf:
            os.utime(tf.name, (0, 0))
            s = [os.path.join(CALC, "v6_bath.py")]
            good = current(tf.name, s, {"v6_bath.py": sha(s[0])})
            stale = current(tf.name, s, {"v6_bath.py": "0"*64}) or current(tf.name, s, {})
        print("selftest: a module imported two levels down is watched:", "yes" if deep else "NO")
        print("selftest: an older log is current only while its source has the recorded hash:",
              "yes" if good and not stale else "NO")
        ok = ok and deep and good and not stale
        sys.exit(0 if ok else 1)
    fails = check(paper, CLAIMS)
    print("v6 numbers: %d claims checked, %d problem(s)" % (len(CLAIMS), len(fails)))
    for f in fails: print("   " + f)
    sys.exit(1 if fails else 0)

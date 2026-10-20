#!/usr/bin/env python3
"""Build INFLATON_v1.md from INFLATON_v1_src.md: number the [@key] citations in order of first
appearance and write the reference list, with the builder of Separate Ways v6 (v6_refs.py).

    python3 pub/inflaton/inflaton_refs.py

Entries shared with v6 come from v6_refs.py. The ones below were read on INSPIRE (titles, authors,
journal, DOI) on 2026-10-08. The build fails if the text cites a key with no entry, or if an entry
below is never cited.
"""
import io, os, re, sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path[:0] = [HERE, os.path.join(HERE, "..", "paper2")]     # in the release v6_refs.py sits beside this file
from v6_refs import REFS as V6, build, CITE

SRC = os.path.join(HERE, "INFLATON_v1_src.md")
OUT = os.path.join(HERE, "INFLATON_v1.md")
EXTRA = {
 "ENT": 'Y. Ema, K. Nakayama, Y. Tang, "Production of purely gravitational dark matter: the case of fermion and vector boson," JHEP 07 (2019) 060, arXiv:1903.10973, https://doi.org/10.1007/JHEP07(2019)060.',
 "Garcia": 'M. A. G. Garcia, Y. Mambrini, K. A. Olive, M. Peloso, "Enhancement of the dark matter abundance before reheating: applications to gravitino dark matter," Phys. Rev. D 96 (2017) 103510, arXiv:1709.01549, https://doi.org/10.1103/PhysRevD.96.103510.',
 "HB": 'N. Herring, D. Boyanovsky, "Gravitational production of nearly thermal fermionic dark matter," Phys. Rev. D 101 (2020) 123522, arXiv:2005.00391, https://doi.org/10.1103/PhysRevD.101.123522.',
 "DFKPP": 'J. F. Dufaux, G. N. Felder, L. Kofman, M. Peloso, D. Podolsky, "Preheating with trilinear interactions: tachyonic resonance," JCAP 07 (2006) 006, hep-ph/0602144, https://doi.org/10.1088/1475-7516/2006/07/006.',
 "EKLRZ": 'K. Enqvist, M. Karčiauskas, O. Lebedev, S. Rusak, M. Zatta, "Postinflationary vacuum instability and Higgs-inflaton couplings," JCAP 11 (2016) 025, arXiv:1608.08848, https://doi.org/10.1088/1475-7516/2016/11/025.',
 "ACT": 'Atacama Cosmology Telescope Collaboration, T. Louis et al., "The Atacama Cosmology Telescope: DR6 power spectra, likelihoods and $\\Lambda$CDM parameters," JCAP 11 (2025) 062, arXiv:2503.14452, https://doi.org/10.1088/1475-7516/2025/11/062.',
 "KL13": 'R. Kallosh, A. Linde, "Universality class in conformal inflation," JCAP 07 (2013) 002, arXiv:1306.5220, https://doi.org/10.1088/1475-7516/2013/07/002.',
 "KLR13": 'R. Kallosh, A. Linde, D. Roest, "Superconformal inflationary $\\alpha$-attractors," JHEP 11 (2013) 198, arXiv:1311.0472, https://doi.org/10.1007/JHEP11(2013)198.',
 "CW73": 'S. R. Coleman, E. J. Weinberg, "Radiative corrections as the origin of spontaneous symmetry breaking," Phys. Rev. D 7 (1973) 1888, https://doi.org/10.1103/PhysRevD.7.1888.',
 "BCGMP": 'N. Bhusal, E. Chávez M., M. A. G. Garcia, A. G. Menkara, M. Pierre, "Fermion (non)reheating with a quartic inflaton potential," JCAP 06 (2026) 064, arXiv:2512.16203, https://doi.org/10.1088/1475-7516/2026/06/064.',
 "KKLSS": 'N. Kaloper, M. Kleban, A. Lawrence, S. Shenker, L. Susskind, "Initial conditions for inflation," JHEP 11 (2002) 037, hep-th/0209231, https://doi.org/10.1088/1126-6708/2002/11/037.',
 "AMM": 'P. R. Anderson, C. Molina-París, E. Mottola, "Short distance and initial state effects in inflation: stress tensor and decoherence," Phys. Rev. D 72 (2005) 043515, hep-th/0504134, https://doi.org/10.1103/PhysRevD.72.043515.',
 "DZ": 'B. Damski, W. H. Zurek, "Adiabatic-impulse approximation for avoided level crossings: from phase transition dynamics to Landau-Zener evolutions and back again," Phys. Rev. A 73 (2006) 063405, cond-mat/0511709, https://doi.org/10.1103/PhysRevA.73.063405.',
 "SAN": 'S. N. Shevchenko, S. Ashhab, F. Nori, "Landau-Zener-Stückelberg interferometry," Phys. Rept. 492 (2010) 1, arXiv:0911.1917, https://doi.org/10.1016/j.physrep.2010.03.002.',
 "VG": 'N. V. Vitanov, B. M. Garraway, "Landau-Zener model: effects of finite coupling duration," Phys. Rev. A 53 (1996) 4288, https://doi.org/10.1103/PhysRevA.53.4288.',
 "GKMO": 'M. A. G. Garcia, K. Kaneta, Y. Mambrini, K. A. Olive, "Reheating and post-inflationary production of dark matter," Phys. Rev. D 101 (2020) 123507, arXiv:2004.08404, https://doi.org/10.1103/PhysRevD.101.123507.',
 "LA": 'K. D. Lozanov, M. A. Amin, "Equation of state and duration to radiation domination after inflation," Phys. Rev. Lett. 119 (2017) 061301, arXiv:1608.01213, https://doi.org/10.1103/PhysRevLett.119.061301.',
 "ENT18": 'Y. Ema, K. Nakayama, Y. Tang, "Production of purely gravitational dark matter," JHEP 09 (2018) 135, arXiv:1804.07471, https://doi.org/10.1007/JHEP09(2018)135.',
 "SW6": 'B. H. Wiseman, "Separate Ways and the Upside Down: corrected graviton-exchange rates and the dark matter of a CPT-symmetric universe," version 6.0, Zenodo (2026), https://doi.org/10.5281/zenodo.22888119.',
 "CGMT": 'M. M. Chaichian, M. Gogberashvili, M. N. Mnatsakanova, T. Tsiskaridze, "CPT violation, mirror world and implications for baryon asymmetry," Eur. Phys. J. C 86 (2026) 425, arXiv:2603.22381, https://doi.org/10.1140/epjc/s10052-026-15592-5.',
 "Kumar": 'N. Kumar, "Creation of a flat universe-antiuniverse pair from nothing," EPL 149 (2025) 69003, https://doi.org/10.1209/0295-5075/adb451.',
 "PS96": 'D. Polarski, A. A. Starobinsky, "Semiclassicality and decoherence of cosmological perturbations," Class. Quant. Grav. 13 (1996) 377, gr-qc/9504030, https://doi.org/10.1088/0264-9381/13/3/006.',
 "KP09": 'C. Kiefer, D. Polarski, "Why do cosmological perturbations look classical to us?," Adv. Sci. Lett. 2 (2009) 164, arXiv:0810.0087, https://doi.org/10.1166/asl.2009.1023.',
 "MD14": 'P. Madau, M. Dickinson, "Cosmic star formation history," Ann. Rev. Astron. Astrophys. 52 (2014) 415, arXiv:1403.0007, https://doi.org/10.1146/annurev-astro-081811-125615.',
 "DESIBAO": 'M. Abdul Karim et al. (DESI Collaboration), "DESI DR2 results. II. Measurements of baryon acoustic oscillations and cosmological constraints," Phys. Rev. D 112 (2025) 083515, arXiv:2503.14738, https://doi.org/10.1103/tr6y-kpc6.',
 "PPlus": 'D. Scolnic et al., "The Pantheon+ analysis: the full data set and light-curve release," Astrophys. J. 938 (2022) 113, arXiv:2112.03863, https://doi.org/10.3847/1538-4357/ac8b7a.',
}

if __name__ == "__main__":
    text = io.open(SRC, encoding="utf-8").read()
    cited = [k for m in CITE.finditer(text) for k in re.findall(r"@([A-Za-z0-9]+)", m.group(1))]
    clash = sorted(set(EXTRA) & set(V6))
    if clash: sys.exit("entries defined twice: %s" % clash)
    pool = dict(V6, **EXTRA)
    refs = {k: pool[k] for k in dict.fromkeys(cited) if k in pool}
    out, missing, _ = build(text, refs)
    unused = [k for k in EXTRA if k not in cited]
    if missing or unused:
        print("missing entries:", missing); print("uncited entries:", unused); sys.exit(1)
    io.open(OUT, "w", encoding="utf-8").write(out)
    print("wrote %s with %d references" % (OUT, len(refs)))

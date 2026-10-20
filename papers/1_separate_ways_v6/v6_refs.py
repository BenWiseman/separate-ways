#!/usr/bin/env python3
"""Build SEPARATE_WAYS_v6.md from SEPARATE_WAYS_v6_src.md: number the [@key] citations in order of
first appearance and write the reference list.

    python3 papers/1_separate_ways_v6/v6_refs.py            # build
    python3 papers/1_separate_ways_v6/v6_refs.py --selftest # the plants must fail

Every entry was read on INSPIRE (titles, authors, journal) on 2026-10-07; equation numbers cited in
the text are those of the arXiv versions named. The build fails if the text cites a key that has no
entry, or if an entry is never cited, so the list cannot drift from the text.
"""
import io, os, re, sys

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "SEPARATE_WAYS_v6_src.md")
OUT = os.path.join(HERE, "SEPARATE_WAYS_v6.md")

REFS = {
 "BFTprl": 'L. Boyle, K. Finn, N. Turok, "CPT-symmetric universe," Phys. Rev. Lett. 121 (2018) 251301, arXiv:1803.08928, https://doi.org/10.1103/PhysRevLett.121.251301.',
 "BFTann": 'L. Boyle, K. Finn, N. Turok, "The Big Bang, CPT, and neutrino dark matter," Annals Phys. 438 (2022) 168767, arXiv:1803.08930, https://doi.org/10.1016/j.aop.2022.168767.',
 "NNP": 'S. Nadal-Gisbert, J. Navarro-Salas, S. Pla, "Low-energy states and CPT invariance at the big bang," Phys. Rev. D 107 (2023) 085018, arXiv:2302.08812, https://doi.org/10.1103/PhysRevD.107.085018.',
 "BTV": 'L. Boyle, N. Turok, V. Vaibhav, "Might the radiation era extend back to the Big Bang? On dark matter production and the relic graviton background in quadratic gravity," arXiv:2609.11750v1, https://doi.org/10.48550/arXiv.2609.11750.',
 "Landau": 'L. D. Landau, "Zur Theorie der Energieübertragung. II," Phys. Z. Sowjetunion 2 (1932) 46.',
 "Zener": 'C. Zener, "Non-adiabatic crossing of energy levels," Proc. R. Soc. Lond. A 137 (1932) 696, https://doi.org/10.1098/rspa.1932.0165.',
 "BNNP": 'P. Beltrán-Palau, S. Nadal-Gisbert, J. Navarro-Salas, S. Pla, "Hadamard and boundary conditions for the Big Bang quantum vacuum," J. Phys. Conf. Ser. 2531 (2023) 012009, arXiv:2204.05414, https://doi.org/10.1088/1742-6596/2531/1/012009.',
 "Fixsen": 'D. J. Fixsen, "The temperature of the cosmic microwave background," Astrophys. J. 707 (2009) 916, arXiv:0911.1955, https://doi.org/10.1088/0004-637X/707/2/916.',
 "Planck": 'Planck Collaboration, N. Aghanim et al., "Planck 2018 results. VI. Cosmological parameters," Astron. Astrophys. 641 (2020) A6, arXiv:1807.06209, https://doi.org/10.1051/0004-6361/201833910.',
 "Garny2016": 'M. Garny, M. Sandora, M. S. Sloth, "Planckian interacting massive particles as dark matter," Phys. Rev. Lett. 116 (2016) 101302, arXiv:1511.03278, https://doi.org/10.1103/PhysRevLett.116.101302.',
 "HLZ": 'T. Han, J. D. Lykken, R.-J. Zhang, "On Kaluza-Klein states from large extra dimensions," Phys. Rev. D 59 (1999) 105006, hep-ph/9811350, https://doi.org/10.1103/PhysRevD.59.105006.',
 "TangWu": 'Y. Tang, Y.-L. Wu, "On thermal gravitational contribution to particle production and dark matter," Phys. Lett. B 774 (2017) 676, arXiv:1708.05138, https://doi.org/10.1016/j.physletb.2017.10.034.',
 "Clery": 'S. Clery, Y. Mambrini, K. A. Olive, S. Verner, "Gravitational portals in the early Universe," Phys. Rev. D 105 (2022) 075005, arXiv:2112.15214, https://doi.org/10.1103/PhysRevD.105.075005.',
 "Bernal": 'N. Bernal, M. Dutra, Y. Mambrini, K. Olive, M. Peloso, M. Pierre, "Spin-2 portal dark matter," Phys. Rev. D 97 (2018) 115020, arXiv:1803.01866, https://doi.org/10.1103/PhysRevD.97.115020.',
 "CEYZ": 'D. J. H. Chung, L. L. Everett, H. Yoo, P. Zhou, "Gravitational fermion production in inflationary cosmology," Phys. Lett. B 712 (2012) 147, arXiv:1109.2524, https://doi.org/10.1016/j.physletb.2012.04.066.',
 "MO": 'Y. Mambrini, K. A. Olive, "Gravitational production of dark matter during reheating," Phys. Rev. D 103 (2021) 115009, arXiv:2102.06214, https://doi.org/10.1103/PhysRevD.103.115009.',
 "DIS": 'S. Davidson, G. Isidori, A. Strumia, "The smallest neutrino mass," Phys. Lett. B 646 (2007) 100, hep-ph/0611389, https://doi.org/10.1016/j.physletb.2007.01.015.',
 "NuFIT": 'I. Esteban, M. C. Gonzalez-Garcia, M. Maltoni, I. Martinez-Soler, J. P. Pinheiro, T. Schwetz, "NuFit-6.0: updated global analysis of three-flavor neutrino oscillations," JHEP 12 (2024) 216, arXiv:2410.05380, https://doi.org/10.1007/JHEP12(2024)216.',
 "DESI": 'W. Elbers et al. (DESI Collaboration), "Constraints on neutrino physics from DESI DR2 BAO and DR1 full shape," Phys. Rev. D 112 (2025) 083513, arXiv:2503.14744, https://doi.org/10.1103/w9pk-xsk7.',
 "BRW": 'C. W. Bauer, N. L. Rodd, B. R. Webber, "Dark matter spectra from the electroweak to the Planck scale," JHEP 06 (2021) 121, arXiv:2007.15001, https://doi.org/10.1007/JHEP06(2021)121.',
 "IceCubeEHE": 'IceCube Collaboration, R. Abbasi et al., "Search for extremely-high-energy neutrinos and first constraints on the ultrahigh-energy cosmic-ray proton fraction with IceCube," Phys. Rev. Lett. 135 (2025) 031001, arXiv:2502.01963, https://doi.org/10.1103/PhysRevLett.135.031001.',
 "Gen2": 'IceCube-Gen2 Collaboration, M. G. Aartsen et al., "IceCube-Gen2: the window to the extreme Universe," J. Phys. G 48 (2021) 060501, arXiv:2008.04323, https://doi.org/10.1088/1361-6471/abbd48.',
 "KM3": 'KM3NeT Collaboration, "Observation of an ultra-high-energy cosmic neutrino with KM3NeT," Nature 638 (2025) 376, https://doi.org/10.1038/s41586-024-08543-1.',
 "AG": 'A. Aguirre, S. Gratton, "Inflation without a beginning: a null boundary proposal," Phys. Rev. D 67 (2003) 083515, gr-qc/0301042, https://doi.org/10.1103/PhysRevD.67.083515.',
 "Neiman": 'Y. Neiman, "Antipodally symmetric gauge fields and higher-spin gravity in de Sitter space," JHEP 10 (2014) 153, arXiv:1406.3291, https://doi.org/10.1007/JHEP10(2014)153.',
 "CG97": 'A. Chamblin, G. W. Gibbons, "Nucleating black holes via nonorientable instantons," Phys. Rev. D 55 (1997) 2177, gr-qc/9607079, https://doi.org/10.1103/PhysRevD.55.2177.',
 "HN": 'D. Harlow, T. Numasawa, "Gauging spacetime inversions in quantum gravity," JHEP 01 (2026) 098, arXiv:2311.09978, https://doi.org/10.1007/JHEP01(2026)098.',
 "companion": 'B. H. Wiseman, "Somewhere Over the Horizon: causal contact under the antipodal identification of horizons, and why it needs four dimensions," preprint, Zenodo (2026), https://doi.org/10.5281/zenodo.23030633.',
 "KT99": 'V. A. Kuzmin, I. I. Tkachev, "Ultrahigh-energy cosmic rays and inflation relics," Phys. Rept. 320 (1999) 199, hep-ph/9903542, https://doi.org/10.1016/S0370-1573(99)00064-2.',
 "CKR": 'D. J. H. Chung, E. W. Kolb, A. Riotto, "Superheavy dark matter," Phys. Rev. D 59 (1998) 023501, hep-ph/9802238, https://doi.org/10.1103/PhysRevD.59.023501.',
 "DMP": 'A. Duran, L. Morrison, S. Profumo, "Sterile neutrino dark matter from generalized CPT-symmetric early-Universe cosmologies," Phys. Rev. D 104 (2021) 023509, arXiv:2103.08626, https://doi.org/10.1103/PhysRevD.104.023509.',
 "Garny2018": 'M. Garny, A. Palessandro, M. Sandora, M. S. Sloth, "Theory and phenomenology of Planckian interacting massive particles as dark matter," JCAP 02 (2018) 027, arXiv:1709.09688, https://doi.org/10.1088/1475-7516/2018/02/027.',
 "Shaposhnikov": 'M. Shaposhnikov, A. Shkerin, I. Timiryasov, S. Zell, "Einstein-Cartan portal to dark matter," Phys. Rev. Lett. 126 (2021) 161301, arXiv:2008.11686, https://doi.org/10.1103/PhysRevLett.126.161301.',
 "ABLMW": 'L. A. Anchordoqui, V. Barger, J. G. Learned, D. Marfatia, T. J. Weiler, "Upgoing ANITA events as evidence of the CPT symmetric universe," LHEP 1 (2018) 13, arXiv:1803.11554, https://doi.org/10.31526/LHEP.1.2018.03.',
 "BDS": 'B. Barman, A. Das, P. Sarmah, "What KM3-230213A event may tell us about the neutrino mass and dark matter," Phys. Rev. D 112 (2025) 075014, arXiv:2504.01447, https://doi.org/10.1103/t58k-sjq5.',
 "Verner": 'S. Verner, "Gravitational production of massive spin-2 particles during reheating," Phys. Rev. D 113 (2026) 123503, arXiv:2512.19041, https://doi.org/10.1103/3st2-jhx5.',
 "KM3decay": 'KM3NeT Collaboration, O. Adriani et al., "Testing heavy dark matter decay as the origin of KM3-230213A," arXiv:2606.09986v1, https://doi.org/10.48550/arXiv.2606.09986.',
 "HV": 'B. L. Hu, E. Verdaguer, "Stochastic gravity: theory and applications," Living Rev. Relativ. 11 (2008) 3, arXiv:0802.0658, https://doi.org/10.12942/lrr-2008-3.',
 "Ishiwata": 'K. Ishiwata, O. Macias, S. Ando, M. Arimoto, "Probing heavy dark matter decays with multi-messenger astrophysical data," JCAP 01 (2020) 003, arXiv:1907.11671, https://doi.org/10.1088/1475-7516/2020/01/003.',
 "PageGeilker": 'D. N. Page, C. D. Geilker, "Indirect evidence for quantum gravity," Phys. Rev. Lett. 47 (1981) 979, https://doi.org/10.1103/PhysRevLett.47.979.',
 "BDOS": 'D. Borah, N. Das, N. Okada, P. Sarmah, "Possible origin of the KM3-230213A neutrino event from dark matter decay," Phys. Rev. D 111 (2025) 123022, arXiv:2503.00097, https://doi.org/10.1103/shsw-mct6.',
 "KolbLong": 'E. W. Kolb, A. J. Long, "Cosmological gravitational particle production and its implications for cosmological relics," Rev. Mod. Phys. 96 (2024) 045005, arXiv:2312.09042, https://doi.org/10.1103/RevModPhys.96.045005.',
 "KT98": 'V. A. Kuzmin, I. I. Tkachev, "Matter creation via vacuum fluctuations in the early universe and observed ultrahigh-energy cosmic ray events," Phys. Rev. D 59 (1999) 123006, hep-ph/9809547, https://doi.org/10.1103/PhysRevD.59.123006.',
 "SWinf": 'B. H. Wiseman, "Upside Down Inflation: dark matter and reheating after a CPT-symmetric bang," preprint, Zenodo (2026).',
 "code": 'B. H. Wiseman, "Separate Ways and the Upside Down: code, figures and the checking pass," version 6.0 [software], Zenodo (2026), https://doi.org/10.5281/zenodo.23034297.',
}

CITE = re.compile(r"\[(@[A-Za-z0-9]+(?:\s*;\s*@[A-Za-z0-9]+)*)\]")

def compress(nums):
    nums = sorted(set(nums)); out = []; i = 0
    while i < len(nums):
        j = i
        while j + 1 < len(nums) and nums[j+1] == nums[j] + 1: j += 1
        out.append(str(nums[i]) if j == i else ("%d,%d" % (nums[i], nums[j]) if j == i + 1 else "%d-%d" % (nums[i], nums[j])))
        i = j + 1
    return ",".join(out)

def build(text, refs):
    order, missing = [], []
    for m in CITE.finditer(text):
        for k in re.findall(r"@([A-Za-z0-9]+)", m.group(1)):
            if k not in refs: missing.append(k)
            elif k not in order: order.append(k)
    unused = [k for k in refs if k not in order]
    num = {k: i + 1 for i, k in enumerate(order)}
    body = CITE.sub(lambda m: "[" + compress([num[k] for k in re.findall(r"@([A-Za-z0-9]+)", m.group(1)) if k in num]) + "]", text)
    lst = "\n\n".join("%d\\. %s" % (num[k], refs[k]) for k in order)
    return body.rstrip() + "\n\n" + lst + "\n", missing, unused

if __name__ == "__main__":
    if "--selftest" in sys.argv:
        _, miss, unused = build("see [@BFTprl; @nosuchkey]", {"BFTprl": "x", "spare": "y"})
        ok = miss == ["nosuchkey"] and unused == ["spare"]
        print("selftest: a missing key and an unused entry are both caught:", "yes" if ok else "NO")
        assert compress([1, 2, 3, 5, 7, 8]) == "1-3,5,7,8"
        sys.exit(0 if ok else 1)
    text = io.open(SRC, encoding="utf-8").read()
    out, missing, unused = build(text, REFS)
    if missing or unused:
        print("missing entries:", missing); print("uncited entries:", unused); sys.exit(1)
    io.open(OUT, "w", encoding="utf-8").write(out)
    print("wrote %s with %d references" % (OUT, len(REFS)))

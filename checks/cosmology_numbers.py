#!/usr/bin/env python3
"""Every high-precision number in the cosmology paper, against every script in both releases.

CLAIMS.tsv anchors the companion, so the cosmology paper's numbers have never had a provenance
gate. The first run of this found two that had drifted from their own script, one of them flagged
inside the script as stale four days earlier, and four in Appendix C.1 with no script at all.

It is NOT in check_all.sh and is not meant to be: it runs 300 scripts and takes the better part of
an hour, against the fifteen minutes the gate suite already costs. Run it before a release.

    python3 checks/cosmology_numbers.py            # report new gaps only
    python3 checks/cosmology_numbers.py --all      # report every unbacked number

The pool is the calculation scripts of both trees and not the figure generators, so a number
stated only in a figure caption can show up unbacked. That is a known edge and not a gap.

THREE THINGS EARLIER VERSIONS GOT WRONG, kept here because all three are easy to repeat.

Scripts must run in their OWN directory. Several in the separate_ways tree source helpers.R by a
relative path, so running them from the repo root produces no output at all and every number they
back looks unbacked. That alone inflated the first report from 18 to 52.

The two trees also want opposite interpreters and directories, and getting it wrong empties the
pool silently rather than loudly. far_side's Python scripts need the repository's .venv, which is
what number_provenance.interpreter() finds when the working directory is the repository root;
called with a bare python3 every one of them dies on numpy, and running from checks/calc hides
the .venv as well. That alone left 32 scripts failing and flagged 39 companion numbers that were
backed all along.

And a bare tolerance is the wrong bar. The paper quotes each figure to its own precision, so a
number agreeing to half a unit in its last quoted place IS reproduced; matching to more digits
than were printed fails on rounding and nothing else.
"""
import io, os, re, sys, glob, subprocess
sys.path.insert(0, os.path.expanduser("~/benlm/tools"))
import number_provenance as np

ROOT  = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PAPER = os.path.join(ROOT, "papers/1_separate_ways/PAPER2_v4_draft.md")

# Known and benign, each with the reason it is not a script's job to produce.
ALLOW = {
    "2511.07517": "an arXiv identifier, not a quantity",
    "2601.19424": "an arXiv identifier, not a quantity",
    "2609.05053": "an arXiv identifier, not a quantity",
    "9.7\\times10^{-48}":  "rho_DM,0, a measured input carried from the cited source",
    "5.966\\times10^{18}": "mu-hat, the Planck-scale mass the production formula carries. NOT the reduced Planck mass, 2.4353e18 GeV: it is sqrt(6) times that, to 1 part in 1e4",
    "0.4142135624": "exactly 1 - sqrt(2); the script checks the agreement, not the digits",
    "4.305": "a log-slope quoted negative; the auditor strips the sign before matching",
    "4.308": "a log-slope quoted negative; the auditor strips the sign before matching",
    "4.006": "a log-slope quoted negative; the auditor strips the sign before matching",
    "1.293\\times10^{-72}": "printed by direct_recoil.R as a comment and a ratio, not a figure",
    "9.588\\times10^{10}":  "the contact-scale floor, stated in D.1 from the displayed formula",
    "1.870\\times10^{13}":  "the bath temperature at H = M_1, from the displayed formula",
}

def collect():
    scripts = sorted(glob.glob(os.path.expanduser("~/separate_ways/**/*.R"), recursive=True)) \
            + sorted(glob.glob(os.path.join(ROOT, "checks/calc/*.R"))) \
            + sorted(glob.glob(os.path.join(ROOT, "checks/calc/*.py")))
    print(f"  running {len(scripts)} scripts...", flush=True)
    outs, failed = [], []
    for sc in scripts:
        # The two trees have opposite conventions and getting this wrong silently empties the
        # pool. separate_ways scripts source helpers.R by a relative path and must run in their
        # own directory. checks/calc scripts are run from the repository root by check_all.sh,
        # which is also where number_provenance looks for the .venv the Python ones need: called
        # with a bare python3 they all die on numpy and every number they back looks unbacked.
        if sc.startswith(ROOT):
            cwd, exe = ROOT, np.interpreter(os.path.relpath(sc, ROOT))
        else:
            cwd = os.path.dirname(sc)
            exe = ["Rscript", os.path.basename(sc)]
        try:
            r = subprocess.run(exe, cwd=cwd, capture_output=True, text=True, timeout=600)
            if r.returncode != 0:
                failed.append(os.path.relpath(sc, os.path.expanduser("~")))
            outs += np.nums_from_out(r.stdout + r.stderr)
        except Exception:
            failed.append(os.path.relpath(sc, os.path.expanduser("~")))
    return outs, failed

CACHE = os.path.join(ROOT, "far_side", ".cosmology_numbers_cache.txt")

def cached_outputs(refresh):
    """Collecting costs the better part of an hour, so it is cached. --refresh recollects."""
    if not refresh and os.path.exists(CACHE):
        vals = [float(x) for x in io.open(CACHE).read().split()]
        print(f"  {len(vals)} numbers from the cache; --refresh to recollect")
        return vals
    outs, failed = collect()
    print(f"  {len(failed)} scripts did not exit clean; {len(outs)} numbers collected")
    io.open(CACHE, "w").write("\n".join(repr(v) for v in outs))
    return outs

def power(outs, hits, dmin, dmax):
    """How often does a WRONG number match anyway? That is one minus the power, and it has to be
    measured rather than assumed: the first version of this check ran at three decimals against a
    pool of twenty-three thousand numbers, where almost any value matches something by chance, so
    its pass carried no information at all. Corrupting each backed number in the band and counting
    how many corruptions still match is the honest statement of what the band is worth."""
    tested = caught = 0
    for lit, _ in hits:
        got = np.nums_from_text(lit)
        if not got: continue
        val, dec, _, sci = got[0]
        if sci or not (dmin <= dec <= dmax): continue
        if not np.matches(val, dec, sci, outs): continue
        tested += 1
        for delta in (7, 3, 9):                       # three corruptions, count it caught if any
            bad = val + delta * 10.0 ** (-dec)
            if not np.matches(bad, dec, sci, outs):
                caught += 1; break
    return tested, caught

def main():
    show_all = "--all" in sys.argv
    body = io.open(PAPER, encoding="utf-8").read().split("\n## References")[0]
    outs = cached_outputs("--refresh" in sys.argv)
    hits = np.uncited(body, min_dec=3)
    new, known = [], 0
    for lit, ctx in hits:
        got = np.nums_from_text(lit)
        if not got:
            continue
        val, dec, _, sci = got[0]
        if np.matches(val, dec, sci, outs):
            continue
        if lit in ALLOW and not show_all:
            known += 1
            continue
        new.append((lit, " ".join(ctx.split())))
    print(f"  {len(hits)} high-precision numbers, {known} unbacked for a recorded reason")

    # What the check is worth, measured. A pass is only evidence where a wrong answer would
    # have failed, and that depends entirely on how many decimals the paper quoted.
    print("     decimals    numbers backed    a corruption is caught    power")
    weak = 0
    for lo, hi, lab in ((3, 3, "3"), (4, 4, "4"), (5, 9, "5+")):
        t, c = power(outs, hits, lo, hi)
        if t == 0:
            print(f"     {lab:>8}    {t:>14}    {'-':>23}    {'-':>5}"); continue
        pw = c / t
        if lo == 3: weak = t
        print(f"     {lab:>8}    {t:>14}    {c:>23}    {pw:>5.0%}")
    print("     A three-decimal match against a pool this size is weak evidence and the band is")
    print("     reported for completeness rather than relied on; the power above says how weak.")

    if not new:
        print("  every other number is reproduced by a script in one of the two releases")
        return 0
    print(f"  {len(new)} NOT reproduced and NOT on the allow list:")
    for lit, ctx in new:
        print(f"     {lit:>20}   ...{ctx[:96]}")
    return 1

if __name__ == "__main__":
    sys.exit(main())

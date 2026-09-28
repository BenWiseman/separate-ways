"""Every term of art gets its plain-English gloss where the reader first meets it.

The pairs live in gloss_terms.tsv. For each one, find the FIRST use of the term in the
manuscript's body and require the gloss to sit within a stated window of it. An edit that
introduces an earlier bare use, or that moves the gloss away from the first use, fails here.
Abstracts are exempt: an arXiv abstract is terse by design and is capped separately.
"""
import re, sys, os

HERE = os.path.dirname(os.path.abspath(__file__))
PUB  = os.path.join(os.path.dirname(HERE), 'pub', 'paper2')

def body_of(name):
    s = open(os.path.join(PUB, name)).read()
    m = re.search(r'^## 1[.\s]', s, flags=re.M)
    a = re.search(r'^## Appendix A', s, flags=re.M)
    if not m:
        raise SystemExit(f'{name}: no "## 1." heading')
    return s[m.start(): a.start() if a else len(s)]

def rows(path):
    out = []
    for line in open(path):
        if not line.strip() or line.lstrip().startswith('#'):
            continue
        f = line.rstrip('\n').split('\t')
        out.append((f[0], f[1], f[2], int(f[3])))
    return out

def check(rs, bodies, quiet=False):
    bad = 0
    for man, term, gloss, win in rs:
        body = bodies[man]
        m = re.search(r'(?<![\w-])' + re.escape(term).replace(r'\ ', r'\s+'), body, flags=re.I)
        if not m:
            print(f'   MISSING  {man:22s} "{term}" does not appear in the body')
            bad += 1
            continue
        lo, hi = max(0, m.start() - 200), m.start() + win
        window = ' '.join(body[lo:hi].split())
        if ' '.join(gloss.split()) in window:
            if not quiet:
                heads = re.findall(r'^#{2,3} .*$', body[:m.start()], flags=re.M)
                head = heads[-1].lstrip('# ')[:42] if heads else '?'
                print(f'   ok       {term:22s} glossed at first use in {head}')
        else:
            print(f'   COLD     {man:22s} "{term}" first used with no gloss within {win} chars')
            print(f'            wanted: {gloss}')
            print(f'            saw:    ...{window[150:400]}...')
            bad += 1
    return bad

def main():
    rs = rows(os.path.join(HERE, 'gloss_terms.tsv'))
    bodies = {n: body_of(n) for n in sorted({r[0] for r in rs})}
    if '--selftest' in sys.argv:
        ok = True
        # a gloss that belongs to a different term must not be found beside this one
        planted = [('PAPER2_v4_draft.md', 'thermofield double',
                    'the two-level problem a solid-state physicist meets', 400)]
        if check(planted, bodies, quiet=True) == 0:
            print('   plant: a term whose gloss is elsewhere is NOT caught -- the check is dead')
            ok = False
        else:
            print('   plant: a term whose gloss sits elsewhere is caught: yes')
        # a term absent from the body must be caught
        if check([('PAPER2_v4_draft.md', 'zzzqqx', 'anything', 10)], bodies, quiet=True) == 0:
            print('   plant: an absent term is NOT caught')
            ok = False
        else:
            print('   plant: an absent term is caught: yes')
        # the real table must pass
        if check(rs, bodies, quiet=True) == 0:
            print('   the real table passes: yes')
        else:
            print('   the real table FAILS')
            ok = False
        return 0 if ok else 1
    bad = check(rs, bodies)
    print(f'   {len(rs)} terms of art, {bad} used cold')
    return 1 if bad else 0

sys.exit(main())

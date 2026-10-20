#!/usr/bin/env python3
"""Separate Ways v6, section 3.1: which orientation class each sheet-exchanging map is in.

A spacetime map is in the CPT class if it reverses time orientation AND total orientation is kept
(in four dimensions: time and all three spatial directions reversed, PT on coordinates), and in the
time-reversal (T) class if it reverses time orientation and total orientation. The CPT theorem
guarantees an antiunitary symmetry only in the first class.
"""
import numpy as np, sys
fails = 0
def check(ok, what):
    global fails
    print("   %-70s %s" % (what, "yes" if ok else "NO  <-- FAILED")); fails += (not ok)
rng = np.random.default_rng(7)

def ds_antipodal_sign(d, trials=200):
    """Sign of the orientation change of X -> -X on dS_d embedded in R^{1,d}."""
    eta = np.diag([-1.0] + [1.0]*d); s = set()
    for _ in range(trials):
        y = rng.normal(size=d); y /= np.linalg.norm(y); t = rng.normal()
        X = np.concatenate([[np.sinh(t)], np.cosh(t)*y])
        _, _, Vt = np.linalg.svd((X @ eta)[None, :]); T = Vt[1:]
        s.add(int(np.sign(np.linalg.det(np.vstack([X, T]))*np.linalg.det(np.vstack([-X, -T])))))
    return s

print("=== the bang map (eta, x) -> (-eta, -x) in four dimensions ===")
B = np.diag([-1.0, -1, -1, -1])
check(np.linalg.det(B) > 0, "keeps total orientation (det = +1): CPT class")
print("\n=== the de Sitter antipodal map X -> -X, by spacetime dimension ===")
for d in [2, 3, 4, 5, 6]:
    s = ds_antipodal_sign(d)
    cls = "keeps orientation (CPT/CRT class)" if s == {1} else "reverses orientation (T class)"
    print("   dS_%d: %s" % (d, cls))
check(ds_antipodal_sign(4) == {-1}, "in four dimensions it reverses orientation: T class")
check(ds_antipodal_sign(5) == {1}, "plant: in five dimensions it keeps orientation, so the test can say yes")
# global coordinates: (tau, n) -> (-tau, -n), and -Id on S^3 in R^4 is a rotation
check(np.linalg.det(-np.eye(4)) > 0, "the antipode of S^3 is a rotation (det -Id_4 = +1), so the map is T up to a rotation")
print("\n=== the antipodal map at a Schwarzschild bifurcation surface ===")
UV = np.diag([-1.0, -1.0]); S2 = None
x = np.array([0, 0, 1.0]); e1 = np.array([1, 0, 0.0]); e2 = np.array([0, 1, 0.0])
s2_sign = np.sign(np.linalg.det([x, e1, e2])*np.linalg.det([-x, -e1, -e2]))
check(np.linalg.det(UV) > 0 and s2_sign < 0, "(U,V) -> -(U,V) keeps the 2D orientation; the S^2 antipode reverses it")
check(np.linalg.det(UV)*s2_sign < 0, "so the Kruskal antipodal map reverses total orientation: T class")
# time orientation: U+V is the timelike direction; it is reversed
Tvec = np.array([0.5, 0.5])            # Kruskal time T = (U+V)/2 as a covector on (U,V)
check(np.allclose(Tvec @ UV, -Tvec), "and it sends Kruskal time T = (U+V)/2 to -T: time orientation reversed")
check(not np.allclose(Tvec @ np.eye(2), -Tvec), "plant: the identity does not reverse T, so that test can fail")
print("\n=== the flat-wedge map J o P_perp ===")
check(np.linalg.det(np.diag([-1.0, -1, 1, 1]) @ np.diag([1.0, 1, -1, -1])) > 0, "diag(-1,-1,1,1) diag(1,1,-1,-1) = -Id_4: CPT class")
sys.exit(1 if fails else 0)

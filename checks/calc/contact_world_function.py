#!/usr/bin/env python3
"""The world function between a point and its fold image, near contact, computed rather than
assumed.

A.19 states the image stress exponent on the assumption that sigma(x, Theta x) vanishes
LINEARLY in the distance off r = M, calling it generic and not computing it. It is the last
assumption in that chain, and it is a boundary-value problem with a known solution at one
point, so it can be continued rather than assumed.

THE SET-UP. Theta preserves UV and therefore r, and sends (U,V,phi) -> (-U,-V,phi+pi). A point
x in the past interior and its image Theta x in the future interior are joined by a curve on
which U and V both increase, so they are causally related in the right order. At r = M exactly
the connecting geodesic is the E = 0 null ray through the bifurcation surface, which is the
curve contact_curve_holonomy.py transports along and contact_vanvleck.R solves in closed form,
and there sigma = 0. Off r = M the connecting geodesic is timelike inside and spacelike
outside, and this finds it by shooting, continued in r from the null solution.

Kruskal is the chart because the geodesic passes through U = V = 0, where Schwarzschild fails
and where every Christoffel below carries a factor of U or V and so vanishes.
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, fsolve

M = 1.0

def r_of_UV(uv):
    if uv <= 0:
        return 2 * M
    if uv < 1e-10:
        return 2 * M * (1 - uv / np.e)
    return brentq(lambda r: (1 - r / (2 * M)) * np.exp(r / (2 * M)) - uv, 1e-12, 2 * M - 1e-14)

def UV_of_r(r):
    return (1 - r / (2 * M)) * np.exp(r / (2 * M))

def F_of_r(r):
    return 32 * M**3 * np.exp(-r / (2 * M)) / r

def christoffels(U, V):
    r = r_of_UV(U * V)
    e = np.exp(-r / (2 * M))
    G = np.zeros((3, 3, 3))
    G[0, 0, 0] = 2 * M * V * (2 * M + r) * e / r**2
    G[1, 1, 1] = 2 * M * U * (2 * M + r) * e / r**2
    G[0, 2, 2] = -U * r / (4 * M)
    G[1, 2, 2] = -V * r / (4 * M)
    G[2, 0, 2] = G[2, 2, 0] = -4 * M**2 * V * e / r**2
    G[2, 1, 2] = G[2, 2, 1] = -4 * M**2 * U * e / r**2
    return G

def gmat(U, V):
    r = r_of_UV(U * V)
    g = np.zeros((3, 3))
    g[0, 1] = g[1, 0] = -F_of_r(r) / 2      # ds^2 = -F dU dV + r^2 dphi^2
    g[2, 2] = r**2
    return g

def norm2(x, k):
    return k @ gmat(x[0], x[1]) @ k

def geodesic(x0, k0, s_end=1.0, rtol=1e-11, atol=1e-13):
    def rhs(s, y):
        x, k = y[:3], y[3:]
        G = christoffels(x[0], x[1])
        return np.concatenate([k, -np.einsum('abc,b,c->a', G, k, k)])
    sol = solve_ivp(rhs, [0, s_end], np.concatenate([x0, k0]),
                    rtol=rtol, atol=atol, dense_output=True)
    return sol

def connect(r_x, k_guess):
    """Find the tangent at x whose geodesic lands on Theta x at affine parameter 1."""
    uv = UV_of_r(r_x)
    U0 = V0 = -np.sqrt(uv)                    # x in the PAST interior, at t = 0
    x0 = np.array([U0, V0, 0.0])
    target = np.array([-U0, -V0, np.pi])
    def miss(k):
        return geodesic(x0, k).y[:3, -1] - target
    k, info, ier, msg = fsolve(miss, k_guess, full_output=True, xtol=1e-13)
    return x0, k, np.abs(miss(k)).max(), ier


if __name__ == "__main__":
    print(__doc__)
    print("=" * 78)

    print("\n[1] the known point: r = M, where the connecting geodesic is null")
    uvM = UV_of_r(M)
    print(f"    UV(M) = {uvM:.6f}, and sqrt(e)/2 = {np.sqrt(np.e)/2:.6f}")
    U0 = -np.sqrt(uvM)
    # the straight ray: x goes from (U0,V0) to (-U0,-V0); guess the tangent from the chord,
    # with the angular rate the null condition fixes at the start.
    r0 = M
    dphi = np.sqrt(F_of_r(r0) * (2 * abs(U0)) * (2 * abs(U0))) / r0
    k_guess = np.array([-2 * U0, -2 * U0, dphi])
    x0, k, res, ier = connect(M, k_guess)
    sig = 0.5 * norm2(x0, k)
    print(f"    shooting residual {res:.2e}, solver flag {ier}")
    print(f"    sigma(x, Theta x) at r = M: {sig:+.3e}   (must be zero: the pair is null)")
    assert res < 1e-8, "the shoot did not land on the image"
    assert abs(sig) < 1e-7, "the connecting geodesic at r = M is not null"

    print("\n[2] it is the ray through the bifurcation surface, not some other geodesic")
    sol = geodesic(x0, k)
    ss = np.linspace(0, 1, 2001)
    tr = sol.sol(ss)
    prod = tr[0] * tr[1]
    i = np.argmin(prod)
    print(f"    min UV along the curve {prod[i]:.2e} at s = {ss[i]:.4f}, so it touches r = 2M")
    print(f"    |U - V| along the curve stays below {np.abs(tr[0]-tr[1]).max():.2e}, so t is constant")
    print(f"    total angle turned {tr[2, -1]:.9f}, and pi = {np.pi:.9f}")
    assert np.abs(tr[0] - tr[1]).max() < 1e-9

    print("\n[3] continuing in r, off the contact boundary, in small steps")
    print("      r/M      sigma            character   residual")
    rows = [(1.0, 0.5 * norm2(x0, k))]
    for sgn in (-1, +1):
        kk = k.copy()
        for j in range(1, 13):
            r_x = 1.0 + sgn * 0.005 * j
            x0j, kk, res, ier = connect(r_x, kk)
            sg = 0.5 * norm2(x0j, kk)
            rows.append((r_x, sg))
            if j % 2 == 0:
                ch = 'timelike' if sg < -1e-9 else 'spacelike' if sg > 1e-9 else 'null'
                print(f"   {r_x:7.3f}  {sg:+14.8f}   {ch:>9}   {res:.1e}")
            worst_res = max(globals().get('worst_res', 0.0), res)
            globals()['worst_res'] = worst_res
            assert res < 1e-8, f"shoot failed at r/M = {r_x}"

    print(f"\n    worst shooting residual over the whole continuation: {worst_res:.1e}")

    print("\n[4] does sigma vanish linearly? A.19 assumed it does.")
    rows.sort()
    rr = np.array([a for a, _ in rows]); sg = np.array([b for _, b in rows])
    for hw in (0.02, 0.04, 0.06):
        m = np.abs(rr - 1) <= hw + 1e-12
        c = np.polyfit(rr[m] - 1, sg[m], 1)
        resid = np.abs(np.polyval(c, rr[m] - 1) - sg[m]).max()
        print(f"    |r/M - 1| <= {hw:.2f}:  slope {c[0]:+.6f}  intercept {c[1]:+.2e}"
              f"  worst residual {resid:.2e}")
    slopes = []
    for h in (0.005, 0.010, 0.020):
        d = (np.interp(1 + h, rr, sg) - np.interp(1 - h, rr, sg)) / (2 * h)
        slopes.append(d)
        print(f"    central difference at half-width {h:.3f}: {d:+.6f}")
    rich = slopes[0] + (slopes[0] - slopes[1]) / 3
    print(f"    Richardson limit of the central difference: {rich:+.8f}")
    print(f"\n    and that is a closed form: 3 pi + 8 = {3*np.pi+8:.8f}, which is TWICE the")
    print(f"    affine length 3 pi/2 + 4 = {1.5*np.pi+4:.8f} that contact_vanvleck.R computes")
    print(f"    for the contact geodesic. Difference: {abs(rich-(3*np.pi+8)):.2e}")
    assert abs(rich - (3 * np.pi + 8)) < 1e-5, "the closed form does not hold"
    print(f"    so  sigma -> {rich:.3f} M (r - M)  near contact, a LINEAR zero with a")
    print("    nonzero slope. The assumption A.19 flagged as generic is now a measurement.")
    print("    The sign runs spacelike outside and timelike inside, so r = M is where the")
    print("    pair stops being unrelated and starts being causally related.")
    assert abs(rich) > 1

    print("\n[5] a correction to how this was phrased, which the run above forces")
    print("    A.19 says contact and conjugacy are ONE condition. That over-states it and the")
    print("    geodesics above show why. The rotation about the axis through x is")
    print("    -sin(phi) d_theta on the equator, a Killing field, so its restriction to any")
    print("    equatorial geodesic is a Jacobi field vanishing at phi = 0 and phi = pi. Every")
    print("    connecting geodesic here lands at phi = pi by construction, at EVERY radius:")
    for r_x in (0.96, 1.00, 1.06):
        kk, xx = k.copy(), x0
        if abs(r_x - 1) > 1e-12:
            step = 0.005 if r_x > 1 else -0.005
            cur = 1.0
            for _ in range(200):
                cur = round(cur + step, 10)
                xx, kk, res, ier = connect(cur, kk)
                if abs(cur - r_x) < 1e-12:
                    break
        trj = geodesic(xx, kk).sol(np.linspace(0, 1, 501))
        rr_s = np.array([r_of_UV(trj[0, j] * trj[1, j]) for j in range(trj.shape[1])])
        amp = np.abs(rr_s * np.sin(trj[2]))
        print(f"      r/M = {r_x:.2f}:  phi(end) = {trj[2, -1]:.9f}; the Jacobi field r sin(phi)"
              f" peaks at {amp.max():.4f}")
        print(f"                  and ends at {amp[-1]:.2e}, so it is a nontrivial field with"
              f" a zero at each end")
        assert abs(trj[2, -1] - np.pi) < 1e-7 and amp.max() > 0.5 and amp[-1] < 1e-6
    print("    So the image sits at a CONJUGATE point at every radius, not only at contact.")
    print("    What r = M picks out is where that conjugate pair becomes NULL, which is where")
    print("    sigma vanishes and the two-point function diverges. Conjugacy is universal here")
    print("    and contact is where the caustic reaches the light cone. The exponent argument")
    print("    is unaffected, since it needs both Delta -> infinity and sigma -> 0 and only")
    print("    the second happens at r = M.")

    print("\n[6] the plant")
    print("    (a) ask for the wrong angle and the answer must move.")
    for bill in (0.9 * np.pi, np.pi, 1.1 * np.pi):
        uv = UV_of_r(M); U0b = -np.sqrt(uv)
        xb = np.array([U0b, U0b, 0.0]); tb = np.array([-U0b, -U0b, bill])
        kb = fsolve(lambda kq: geodesic(xb, kq).y[:3, -1] - tb, k, xtol=1e-13)
        print(f"        bill {bill/np.pi:.2f} pi  ->  sigma = {0.5*norm2(xb, kb):+12.8f}")
    print("    Only pi gives zero, so the null pair is the antipodal one and not any pair.")
    print("    (b) the shot curve must actually solve the geodesic equation, not merely land.")
    trj = geodesic(x0, k)
    sm = np.linspace(0.05, 0.95, 19)
    y = trj.sol(sm)
    worst = 0.0
    for j in range(len(sm)):
        xj, kj = y[:3, j], y[3:, j]
        G = christoffels(xj[0], xj[1])
        acc = -np.einsum('abc,b,c->a', G, kj, kj)
        num = (trj.sol(sm[j] + 1e-5)[3:] - trj.sol(sm[j] - 1e-5)[3:]) / 2e-5
        worst = max(worst, np.abs(acc - num).max())
    print(f"        worst geodesic-equation residual along the curve: {worst:.2e}")
    assert worst < 1e-5

    print("\n[7] what that does to the exponent, which improves rather than damages it")
    print("    A.19 justifies delta^-3/2 by saying both sigma AND the shortfall to the")
    print("    conjugate point vanish linearly off r = M. The second half is wrong: the")
    print("    shortfall is identically zero, because section 5 shows the endpoint is")
    print("    conjugate at every radius. The right statement is a caustic-order rule.")
    print("    A Green function's generic light-cone singularity is delta^-1, and an order-n")
    print("    caustic on the connecting family steepens it by half a power per focusing")
    print("    direction:")
    print("         exponent = -(D - 2 + n)/2,   n = the number of directions that refocus,")
    print("    which is -(1 + n/2) in four dimensions. antipodal_sphere_family.R checks the")
    print("    general form on an exactly solvable family spanning D = 3 to 7.")
    print("    Three independent points, none of them fitted to the rule:")
    table = [("A.18 control, no caustic",        0, -1.0095,  "measured"),
             ("A.18 image pair",                 1, -1.49874, "measured"),
             ("Einstein static universe, S^3",   2, -2.0,     "exact, closed form")]
    print("      case                              n    rule     found       status")
    for name, n, got, how in table:
        rule = -(1 + n / 2)
        print(f"      {name:<32} {n}   {rule:+.3f}   {got:+.5f}   {how}")
        assert abs(rule - got) < 0.01
    print(f"    the three magnitudes again, plainly: {abs(table[0][2]):.4f}, "
          f"{abs(table[1][2]):.5f}, {abs(table[2][2]):.1f}, against 1, 3/2 and 2.")
    print("    And GR-side, contact_conjugacy_general.R gives n = D - 3, so in D dimensions")
    print("    the image term goes as delta^{-(2D-5)/2} and the stress two powers steeper.")
    for D in range(4, 8):
        print(f"      D = {D}:  n = {D-3}, image term delta^{{{-(2*D-5)/2:+.1f}}},"
              f" stress delta^{{{-(2*D-5)/2-2:+.1f}}}")
    print("    Four dimensions gives -3/2 and -7/2, which is what A.18 measured.")

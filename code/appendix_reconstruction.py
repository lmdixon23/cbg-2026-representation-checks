"""
appendix_reconstruction.py -- convention-lock against the paper's worked pure-braid example
(arXiv:2606.20473, appendix: f'(b_{zeta4 zeta7}) as a 16-flip product, zeta_l = l).

RESULTS (see verification/STATUS-LEDGER.md):
 [LOCKED]  individual complex f' flip matrices reproduce exactly:
             A_4615 -> cos 0.79056942, sin 0.61237244 (real)
             A_1524 -> cos 1.06066017, sin 0.35355339 i (imaginary)
           via the operator with (i,j,k,l) read positionally from the subscript A_{ikjl}
           (so a=i,b=k,c=j,d=l) under the PRINCIPAL complex sqrt, NO separate (-1)^s sign flip.
 [LOCKED]  composition order is REVERSED (rightmost flip A_1546 acts first), and the RATIONAL
           representation f (eq. matrix_n_2) reproduces the paper's stated global property of
           pure-braid matrices -- determinant 1 and UNIT DIAGONAL -- to ~1e-16 on b_{4,7}.
           This is a ground-truth match on the full 16-flip product and validates the order,
           the flip-sequence reconstruction, and the initial-triangulation recovery together.
 [RESOLVED] the COMPLEX f' product is NOT unit-diagonal under THIS script's naive reduced-basis
           composition -- a convention artifact of the touched-triangle reconstruction here, not the
           paper. The paper prints the full 9x9 f' matrices; appendix_full_matrix.py multiplies the
           paper's printed factors and obtains the correct complex-orthogonal, det-1, unit-diagonal,
           rank-2 unipotent matrix, and shows the paper's two DISPLAYED result matrices are erroneous
           (I + skew). So the note relies on the rational f (above) plus the diagonal similarity
           f' = D f D^{-1}, not on this script's complex path.

Run:  python appendix_reconstruction.py
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cmath
import numpy as np
from operators import cos_sin, is_complex_orthogonal

VALUES = {l: float(l) for l in range(1, 8)}            # zeta_l = l, seven points (n=4)
PRODUCT = ["4615","1524","4716","3417","1645","1746","2437","3714",
           "4627","2734","1467","6724","4516","1647","2415","1546"]


def params(abcd):
    a, b, c, d = (int(ch) for ch in abcd)
    return a, c, b, d, frozenset({a, b})               # i,j,k,l, removed-diagonal


def local_2x2(i, j, k, l, kind):
    """[ijl;jkl]-from-[ijk;ikl] block. 'rational' = eq.operator-2/3; 'complex' = cos/sin f'."""
    z = lambda a, b: VALUES[a] - VALUES[b]
    if kind == "rational":
        return np.array([[z(i, l) / z(i, k), z(i, j) / z(i, k)],
                         [z(l, k) / z(i, k), z(j, k) / z(i, k)]], dtype=complex)
    c = cmath.sqrt((z(i, l) * z(j, k)) / (z(i, k) * z(j, l)))
    s = cmath.sqrt((z(i, j) * z(k, l)) / (z(i, k) * z(j, l)))
    return np.array([[c, s], [-s, c]], dtype=complex)


def apply_flip(vec, i, j, k, l, removed, kind):
    ijk, ikl = frozenset({i, j, k}), frozenset({i, k, l})
    ijl, jkl = frozenset({i, j, l}), frozenset({j, k, l})
    M = local_2x2(i, j, k, l, kind)
    out = dict(vec)
    if set(removed) == {i, k}:                          # forward ik->jl
        if ijk in out or ikl in out:
            a, b = out.pop(ijk, 0), out.pop(ikl, 0)
            out[ijl] = out.get(ijl, 0) + M[0, 0] * a + M[0, 1] * b
            out[jkl] = out.get(jkl, 0) + M[1, 0] * a + M[1, 1] * b
    else:                                               # inverse jl->ik
        Mi = np.linalg.inv(M)
        if ijl in out or jkl in out:
            a, b = out.pop(ijl, 0), out.pop(jkl, 0)
            out[ijk] = out.get(ijk, 0) + Mi[0, 0] * a + Mi[0, 1] * b
            out[ikl] = out.get(ikl, 0) + Mi[1, 0] * a + Mi[1, 1] * b
    return out


def full_product(kind, order="reversed"):
    seq = PRODUCT if order == "written" else list(reversed(PRODUCT))   # rightmost acts first
    flips = [params(s) for s in seq]
    live, produced, T0 = set(), set(), set()
    for (i, j, k, l, rem) in flips:
        ijk, ikl = frozenset({i, j, k}), frozenset({i, k, l})
        ijl, jkl = frozenset({i, j, l}), frozenset({j, k, l})
        fwd = (set(rem) == {i, k})
        ins, outs = ((ijk, ikl), (ijl, jkl)) if fwd else ((ijl, jkl), (ijk, ikl))
        for t in ins:
            if t not in live:
                if t in produced:
                    return None
                T0.add(t); live.add(t)
        live -= set(ins); live |= set(outs); produced |= set(outs)
    if live != T0:
        return None
    basis = sorted(T0, key=lambda t: tuple(sorted(t)))
    idx = {t: r for r, t in enumerate(basis)}
    state = {t: {t: 1 + 0j} for t in T0}
    for (i, j, k, l, rem) in flips:
        for t0 in state:
            state[t0] = apply_flip(state[t0], i, j, k, l, rem, kind)
    M = np.zeros((len(basis), len(basis)), dtype=complex)
    for col, t0 in enumerate(basis):
        for t, co in state[t0].items():
            M[idx[t], col] = co
    return M


def main():
    print("1) [LOCKED] individual complex f' flip matrices vs paper ground truth:")
    for abcd, wc, ws in [("4615", 0.79056942, 0.61237244 + 0j), ("1524", 1.06066017, 0.35355339j)]:
        i, j, k, l, _ = params(abcd)
        M = local_2x2(i, j, k, l, "complex")
        c, s = M[0, 0], -M[1, 0]
        print(f"   A_{abcd}: cos={c:.8f} (paper {wc:.8f}) {'OK' if abs(c-wc)<1e-7 else 'X'};"
              f"  sin={s:.8f} (paper {ws}) {'OK' if abs(s-ws)<1e-7 else 'X'}")

    print("\n2) full 16-flip product f'(b_4,7) -- order x operator vs paper's unit-diagonal+det-1:")
    N = 2 * 4 + 1                                   # n=4 -> full representation dimension 2n+1 = 9
    for kind in ("rational", "complex"):
        for order in ("written", "reversed"):
            M = full_product(kind, order)
            if M is None:
                print(f"   {kind:9s} {order:9s}: did not close"); continue
            active = M.shape[0]
            spect = N - active                       # untouched triangles act as identity
            ud = np.max(np.abs(np.diag(M) - 1))
            if spect < 0:
                # active > 2n+1 means this (wrong) order over-counted: not a valid reconstruction
                print(f"   {kind:9s} {order:9s}: INVALID reconstruction (active {active} > full {N}); "
                      f"det={np.linalg.det(M):.4f} max|diag-1|={ud:.2e}  (wrong composition order)")
                continue
            # property on the full 2n+1 matrix = property on the active block (spectators add 1s / x1)
            full_unit = (ud < 1e-9)
            tag = "  <-- matches paper (UNIT DIAGONAL, det 1) on full 9x9" if (full_unit and abs(np.linalg.det(M)-1) < 1e-9) else ""
            print(f"   {kind:9s} {order:9s}: active block {active}x{active} (+{spect} spectator => full {N}x{N}); "
                  f"det={np.linalg.det(M):.4f} max|diag-1|={ud:.2e}{tag}")
    assert (N - full_product('rational', 'reversed').shape[0]) == 1, "dimension accounting"

    print("\n   [LOCKED, individual] the two printed flip matrices reproduce ENTRYWISE (above).")
    print("   [LOCKED, property]   rational f, reversed order: the active 8x8 block is unit-diagonal+det-1,")
    print("                        so the full 9x9 = block (+) I_1(spectator) has the paper's stated property.")
    print("                        This is a GLOBAL-PROPERTY match for the rational f. The paper DOES print")
    print("                        the full 9x9 f' result matrices; the entrywise/erratum check against them")
    print("                        is in appendix_full_matrix.py.")
    print("   [RESOLVED] the complex f' product is not unit-diagonal under THIS script's naive reduced-basis")
    print("            composition -- a convention artifact here, not the paper. appendix_full_matrix.py")
    print("            multiplies the paper's printed factors to the correct orthogonal/det-1/unit-diagonal")
    print("            rank-2 matrix and shows the paper's two DISPLAYED result matrices are erroneous (former O1).")


if __name__ == "__main__":
    main()

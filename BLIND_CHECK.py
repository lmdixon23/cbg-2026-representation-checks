"""
BLIND_CHECK.py -- the author's independent recomputation.

Written WITHOUT importing anything from code/. Everything below is recomputed from scratch from the
operator formulas of arXiv:2606.20473 and zeta_l = l, then compared against the paper's printed values
and the note's stated properties. If this script and code/ ever disagree, stop and report the
disagreement; do not edit either to match the other.

Recomputes, independently:
  (1) the two printed appendix flip matrices A_4615, A_1524 (paper ground truth);
  (2) the Ptolemy/Plucker identity and hence cos^2 + sin^2 = 1, over all orderings;
  (3) the single rational flip preserves the diagonal form Q (f^T Q_p f = Q_t), symbolically;
  (4) the worked pure braid b_4,7 (rational f, 16 flips, reversed order): determinant 1, unit
      diagonal, rank-2 one-step unipotent.

Run:  python BLIND_CHECK.py     (a few seconds)
Exit code 0 and "BLIND CHECK PASSED" on success; an AssertionError otherwise.
"""
import cmath
import itertools
import numpy as np
import sympy as sp

# ---------------------------------------------------------------- zeta_l = l, seven points (n=4)
VALUES = {l: float(l) for l in range(1, 8)}
def z(a, b):
    return VALUES[a] - VALUES[b]

def cos_sin(i, j, k, l):
    """Paper's cos/sin for the flip ik->jl (principal complex sqrt)."""
    c = cmath.sqrt((z(i, l) * z(j, k)) / (z(i, k) * z(j, l)))
    s = cmath.sqrt((z(i, j) * z(k, l)) / (z(i, k) * z(j, l)))
    return c, s

# Subscript A_{abcd} -> operator indices read positionally: a=i, b=k, c=j, d=l; diagonal {a,b}={i,k}.
def indices(abcd):
    a, b, c, d = (int(ch) for ch in abcd)
    return a, c, b, d            # (i, j, k, l)


# ---------------------------------------------------------------- (1) the two printed flip matrices
def check_printed_flip_matrices():
    ground_truth = {"4615": (0.79056942, 0.61237244 + 0j),
                    "1524": (1.06066017, 0.35355339j)}
    print("(1) printed appendix flip matrices vs paper ground truth:")
    ok = True
    for abcd, (wc, ws) in ground_truth.items():
        i, j, k, l = indices(abcd)
        c, s = cos_sin(i, j, k, l)
        mc, ms = abs(c - wc) < 1e-7, abs(s - ws) < 1e-7
        ok = ok and mc and ms
        print(f"    A_{abcd} (i,j,k,l)={i,j,k,l}: cos={c:.8f} (paper {wc}) {'OK' if mc else 'X'}; "
              f"sin={s:.8f} (paper {ws}) {'OK' if ms else 'X'}")
    assert ok, "a printed flip matrix did not reproduce from the formulas"
    return ok


# ---------------------------------------------------------------- (2) Ptolemy => cos^2 + sin^2 = 1
def check_ptolemy_orthogonality():
    worst = 0.0
    for perm in itertools.permutations([2.0, 3.0, 5.0, 8.0]):
        global VALUES
        save = VALUES
        VALUES = {1: perm[0], 2: perm[1], 3: perm[2], 4: perm[3]}
        c, s = cos_sin(1, 2, 3, 4)
        worst = max(worst, abs(c * c + s * s - 1))
        VALUES = save
    print(f"(2) cos^2 + sin^2 = 1 (Ptolemy) over all 24 orderings: max |residual| = {worst:.2e}")
    assert worst < 1e-12, "Ptolemy/orthogonality failed"
    return True


# ---------------------------------------------------------------- (3) rational flip preserves Q (exact)
def check_Q_isometry_symbolic():
    zi, zj, zk, zl = sp.symbols("zi zj zk zl")
    f = sp.Matrix([[(zi - zl) / (zi - zk), (zi - zj) / (zi - zk)],
                   [(zl - zk) / (zi - zk), (zj - zk) / (zi - zk)]])   # cols D_ijk,D_ikl -> rows D_ijl,D_jkl
    Q = lambda a, b, c: 1 / ((a - b) * (a - c) * (b - c))
    Qt = sp.diag(Q(zi, zj, zk), Q(zi, zk, zl))                       # source (D_ijk, D_ikl)
    Qp = sp.diag(Q(zi, zj, zl), Q(zj, zk, zl))                       # target (D_ijl, D_jkl)
    residual = sp.simplify(f.T * Qp * f - Qt)
    ok = residual == sp.zeros(2, 2)
    print(f"(3) single rational flip preserves Q=1/(z_pq z_pr z_qr):  f^T Q_p f - Q_t = 0 (exact): {ok}")
    assert ok, "single-flip Q-isometry identity failed"
    return ok


# ---------------------------------------------------------------- (4) b_4,7 rational: det 1, unit diag
PRODUCT_B47 = ["4615", "1524", "4716", "3417", "1645", "1746", "2437", "3714",
               "4627", "2734", "1467", "6724", "4516", "1647", "2415", "1546"]

def b47_rational_matrix():
    """Compose the 16 rational flips (rightmost acts first) on the touched-triangle basis.
    Every flip here removes the diagonal {a,b}={i,k}, so all flips are 'forward'."""
    flips = [indices(s) for s in reversed(PRODUCT_B47)]
    tri = lambda *v: frozenset(v)
    live, produced, T0 = set(), set(), set()
    for (i, j, k, l) in flips:
        ins = (tri(i, j, k), tri(i, k, l))
        outs = (tri(i, j, l), tri(j, k, l))
        for t in ins:
            if t not in live:
                T0.add(t); live.add(t)
        live -= set(ins); live |= set(outs); produced |= set(outs)
    basis = sorted(T0, key=lambda t: tuple(sorted(t)))
    idx = {t: r for r, t in enumerate(basis)}
    state = {t: {t: 1 + 0j} for t in T0}                            # image of each initial basis triangle
    for (i, j, k, l) in flips:
        ijk, ikl = tri(i, j, k), tri(i, k, l)
        ijl, jkl = tri(i, j, l), tri(j, k, l)
        A, B = z(i, l) / z(i, k), z(i, j) / z(i, k)
        C, D = z(l, k) / z(i, k), z(j, k) / z(i, k)
        for v in state.values():
            a = v.pop(ijk, 0); b = v.pop(ikl, 0)
            if a or b:
                v[ijl] = v.get(ijl, 0) + A * a + B * b
                v[jkl] = v.get(jkl, 0) + C * a + D * b
    M = np.zeros((len(basis), len(basis)), dtype=complex)
    for col, t0 in enumerate(basis):
        for t, co in state[t0].items():
            M[idx[t], col] = co
    return M

def check_b47_property():
    M = b47_rational_matrix()
    N = M - np.eye(M.shape[0])
    det = np.linalg.det(M)
    unit = np.max(np.abs(np.diag(M) - 1))
    rank = np.linalg.matrix_rank(N, tol=1e-9)
    nilp = np.allclose(N @ N, 0, atol=1e-9)
    print(f"(4) b_4,7 rational f ({M.shape[0]}x{M.shape[0]} active block): det={det.real:.6f}, "
          f"max|diag-1|={unit:.2e}, rank(M-I)={rank}, (M-I)^2=0={nilp}")
    assert abs(det - 1) < 1e-9 and unit < 1e-9 and rank == 2 and nilp, \
        "b_4,7 rational matrix is not the expected det-1 unit-diagonal rank-2 unipotent"
    return True


def main():
    print("== BLIND CHECK (independent recomputation; no code/ imports) ==\n")
    check_printed_flip_matrices()
    check_ptolemy_orthogonality()
    check_Q_isometry_symbolic()
    check_b47_property()
    print("\nBLIND CHECK PASSED -- printed flip matrices, Ptolemy/orthogonality, the Q-isometry")
    print("identity, and the b_4,7 determinant/unit-diagonal/rank-2 properties all reproduce.")


if __name__ == "__main__":
    main()

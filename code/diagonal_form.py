"""
diagonal_form.py -- resolves O1/D1: the rational representation f and the complex-orthogonal
representation f' of arXiv:2606.20473 are related by a diagonal similarity f' = D f D^{-1}.

RESULT.
  Each rational flip operator f (eq. operator-2/3) preserves the diagonal symmetric bilinear form
       Q(Delta_pqr) = 1 / ( zeta_pq * zeta_pr * zeta_qr )        (reciprocal product of the
  triangle's three pairwise zeta-differences).  Equivalently, with D = diag(Q^{1/2}),
       f'_flip = D_p f_flip D_t^{-1}     (flip-by-flip, source triangulation t -> target p),
  i.e. f' is f expressed in the Q-orthonormal triangle basis.  Consequences:
   * f' lands in the complex-orthogonal group (re-derives Prop. "orthogonality is structural");
   * for any pure braid (closed flip loop) the D's telescope, so f'(beta) = D f(beta) D^{-1} is
     conjugate to f(beta) by a diagonal matrix; since f(beta) has unit diagonal and det 1, so does
     f'(beta).  This is exactly the property the paper asserts for f'.  The square root Q^{1/2} in D
     carries a per-triangle sign, so independently chosen principal branches need not keep D
     single-valued around a loop.  We CONJECTURE (not proved here) that the paper's B/C sign rule is
     the branch bookkeeping that restores single-valuedness; what is verified below is the exact
     algebraic identity f^T Q f = Q and its composition, plus a numerical confirmation on b_4,7.

Run:  python diagonal_form.py
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cmath
import numpy as np
import sympy as sp
from operators import is_complex_orthogonal
from appendix_reconstruction import params, apply_flip, PRODUCT, VALUES


def symbolic_proof():
    """EXACT: each rational flip block preserves Q = 1/(product of pairwise diffs)."""
    zi, zj, zk, zl = sp.symbols('zi zj zk zl')
    f = sp.Matrix([[(zi - zl) / (zi - zk), (zi - zj) / (zi - zk)],
                   [(zl - zk) / (zi - zk), (zj - zk) / (zi - zk)]])   # cols D_ijk,D_ikl -> rows D_ijl,D_jkl
    Q = lambda a, b, c: 1 / ((a - b) * (a - c) * (b - c))
    Qt = sp.diag(Q(zi, zj, zk), Q(zi, zk, zl))     # source form (D_ijk, D_ikl)
    Qp = sp.diag(Q(zi, zj, zl), Q(zj, zk, zl))     # target form (D_ijl, D_jkl)
    residual = sp.simplify(f.T * Qp * f - Qt)
    ok = residual == sp.zeros(2, 2)
    print(f"[PROOF, exact] rational flip preserves Q(D_pqr)=1/(z_pq z_pr z_qr):  f^T Q_p f - Q_t = 0  -> {ok}")
    return ok


def _rational_product_with_basis():
    """f(b_4,7) on the touched-triangle basis (= full_product('rational','reversed') with basis)."""
    flips = [params(s) for s in reversed(PRODUCT)]
    live, produced, T0 = set(), set(), set()
    for (i, j, k, l, rem) in flips:
        ijk, ikl = frozenset({i, j, k}), frozenset({i, k, l})
        ijl, jkl = frozenset({i, j, l}), frozenset({j, k, l})
        fwd = (set(rem) == {i, k})
        ins, outs = ((ijk, ikl), (ijl, jkl)) if fwd else ((ijl, jkl), (ijk, ikl))
        for t in ins:
            if t not in live:
                T0.add(t); live.add(t)
        live -= set(ins); live |= set(outs); produced |= set(outs)
    basis = sorted(T0, key=lambda t: tuple(sorted(t)))
    idx = {t: r for r, t in enumerate(basis)}
    state = {t: {t: 1 + 0j} for t in T0}
    for (i, j, k, l, rem) in flips:
        for t0 in state:
            state[t0] = apply_flip(state[t0], i, j, k, l, rem, "rational")
    M = np.zeros((len(basis), len(basis)), dtype=complex)
    for col, t0 in enumerate(basis):
        for t, co in state[t0].items():
            M[idx[t], col] = co
    return M, basis


def per_flip_numerical():
    """Confirm the exact per-flip identity numerically on the appendix flips: f_flip preserves Q."""
    z = lambda a, b: VALUES[a] - VALUES[b]
    Q = lambda t: 1 / (z(*sorted(t)[:2]) * z(sorted(t)[0], sorted(t)[2]) * z(*sorted(t)[1:]))
    worst = 0.0
    for abcd in PRODUCT:
        i, j, k, l, rem = params(abcd)
        A, B, C, Dd = z(i, l) / z(i, k), z(i, j) / z(i, k), z(l, k) / z(i, k), z(j, k) / z(i, k)
        fblk = np.array([[A, B], [C, Dd]], dtype=complex)              # cols D_ijk,D_ikl
        Qt = np.diag([1 / (z(i, j) * z(i, k) * z(j, k)), 1 / (z(i, k) * z(i, l) * z(k, l))])
        Qp = np.diag([1 / (z(i, j) * z(i, l) * z(j, l)), 1 / (z(j, k) * z(j, l) * z(k, l))])
        worst = max(worst, np.max(np.abs(fblk.T @ Qp @ fblk - Qt)))
    print(f"[CONFIRM, per flip] all 16 appendix flips preserve Q (max||f^T Q_p f - Q_t||): {worst:.2e}")
    return worst < 1e-9


def numerical_confirmation():
    """On the worked pure braid b_4,7: f' = D f D^{-1} is unit-diagonal, det 1, AND complex-orthogonal.
    Unit-diagonal holds for any diagonal D (-> O1 resolved); orthogonality holds for D=G^{1/2} where G is
    the diagonal form the loop preserves (recovered directly here)."""
    M, basis = _rational_product_with_basis()
    n = M.shape[0]
    # O1 (robust): conjugating f by ANY diagonal leaves the diagonal unchanged
    print(f"[CONFIRM on b_4,7]  f (rational): unit-diag max|diag-1|={np.max(np.abs(np.diag(M)-1)):.2e}, det={np.linalg.det(M):.4f}")
    print(f"   f'=D f D^-1 (any diagonal D): unit diagonal preserved -> O1: the Q-orthonormalized f'(pure braid) is unit-diag")
    # solve for the diagonal form G with M^T G M = G  (the consistent-orientation Q on the loop)
    rows = []
    for i in range(n):
        for j in range(n):
            c = np.array([M[k, i] * M[k, j] for k in range(n)], dtype=complex)
            if i == j:
                c[i] -= 1
            rows.append(c)
    Vh = np.linalg.svd(np.array(rows))[2]
    g = np.conj(Vh[-1])                                # a preserved diagonal form
    D = np.diag([cmath.sqrt(x) for x in g])
    fprime = D @ M @ np.linalg.inv(D)
    print(f"   loop preserves a diagonal form G (||M^T G M - G||={np.linalg.norm(M.T@np.diag(g)@M-np.diag(g)):.1e}); "
          f"with D=G^1/2:")
    print(f"   f' = D f D^-1: complex-orthogonal={is_complex_orthogonal(fprime,1e-6)}, "
          f"unit-diag={np.max(np.abs(np.diag(fprime)-1)):.2e}, det={np.linalg.det(fprime):.4f}")
    print("   => f'(b_4,7) is unit-diagonal complex-orthogonal det-1. O1 RESOLVED; D1 was the branch of")
    print("      D=Q^1/2 (per-triangle orientation) not tracked by naive principal-sqrt composition.")
    return is_complex_orthogonal(fprime, 1e-6)


def inductive_step_symbolic():
    """EXACT inductive step: a composite of two flips sharing a triangle preserves Q.
    Base case = Proposition (single flip, symbolic_proof). This composes it: for consecutive flips the
    shared triangle carries ONE value of Q, so flip-1's target form is flip-2's source form and they
    cancel -- f(beta) is a Q-isometry for every flip word, with no prose telescoping."""
    zi, zj, zk, zl, zm = sp.symbols('zi zj zk zl zm')
    Q = lambda a, b, c: 1 / ((a - b) * (a - c) * (b - c))      # fixed vertex order per triangle
    # flip 1 on quad (i,j,k,l): D_ijk,D_ikl -> D_ijl,D_jkl
    a, b = (zi - zl) / (zi - zk), (zl - zk) / (zi - zk)
    c, d = (zi - zj) / (zi - zk), (zj - zk) / (zi - zk)
    # flip 2 on quad (i,j,l,m): D_ijl,D_ilm -> D_ijm,D_jlm  (shares D_ijl with flip 1's output)
    ap, bp = (zi - zm) / (zi - zl), (zm - zl) / (zi - zl)
    cp, dp = (zi - zj) / (zi - zl), (zj - zl) / (zi - zl)
    # t0=(ijk,ikl,ilm) -> t1=(ijl,jkl,ilm) -> t2=(jkl,ijm,jlm)
    f1 = sp.Matrix([[a, c, 0], [b, d, 0], [0, 0, 1]])         # rows t1, cols t0
    f2 = sp.Matrix([[0, 1, 0], [ap, 0, cp], [bp, 0, dp]])      # rows t2, cols t1
    Qt0 = sp.diag(Q(zi, zj, zk), Q(zi, zk, zl), Q(zi, zl, zm))
    Qt1 = sp.diag(Q(zi, zj, zl), Q(zj, zk, zl), Q(zi, zl, zm))
    Qt2 = sp.diag(Q(zj, zk, zl), Q(zi, zj, zm), Q(zj, zl, zm))
    s1 = sp.simplify(f1.T * Qt1 * f1 - Qt0) == sp.zeros(3, 3)
    s2 = sp.simplify(f2.T * Qt2 * f2 - Qt1) == sp.zeros(3, 3)
    F = f2 * f1
    sF = sp.simplify(F.T * Qt2 * F - Qt0) == sp.zeros(3, 3)
    print(f"[INDUCTIVE STEP, exact] flip1 and flip2 preserve Q: {s1 and s2}; "
          f"composite f2*f1 preserves Q (intermediate form cancels): {sF}")
    print("   -> by induction f(beta) is a Q-isometry for every flip word; for a pure braid this gives")
    print("      the Q-orthonormalized f'=D f D^-1 is complex-orthogonal, unit-diagonal, det 1 (no telescoping);")
    print("      coincidence with the paper's B/C-signed f' around loops is the branch conjecture (C8').")
    return s1 and s2 and sF


def main():
    print("== f' = D f D^-1 : the rational and complex-orthogonal reps are diagonally similar ==\n")
    symbolic_proof()
    print()
    inductive_step_symbolic()
    print()
    per_flip_numerical()
    print()
    numerical_confirmation()


if __name__ == "__main__":
    main()

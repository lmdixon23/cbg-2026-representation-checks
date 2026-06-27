"""
exact_kernel.py -- EXACT (sympy, no floating point) lift of the C9 abelianization structure for the
inner pure braid group PB_n, n = 3, 4.  Promotes float_kernel.py's numerical finding to exact arithmetic.

Method: the braid word's *combinatorial* flip sequence is detected from the planar geometry (float is
used only to identify which discrete 2->2 flips occur, in which order), then the EXACT rational flip
operators (eq. operator-2/3, with zeta_l = l so all entries are rationals) are applied symbolically.
By C8 the per-flip diagonal similarity f'=DfD^{-1} holds (diagonal_form.py); transfer of this structure to
the paper's B/C-signed f' is conditional on the C8' branch conjecture, so the EXACT result below is for the
rational f (ker f' is OPEN).

Verified EXACTLY for n = 3 and n = 4:
  * each inner generator A_ij -> I + N_ij with N_ij^2 = 0 and rank N_ij = 2 (rank-2 one-step unipotent)
  * N_ij N_kl = 0 for ALL pairs (products vanish, not merely brackets) -> the image is abelian
  * the N_ij are linearly independent, dim = C(n,2) = dim PB_n^{ab}
  => f restricted to inner PB_n is the abelianization; ker f = [PB_n, PB_n], in EXACT arithmetic GIVEN the
     Delaunay-located flip sequences (ker f' is OPEN, conditional on the C8' branch conjecture).

Run:  python exact_kernel.py   (n=4 takes ~1-2 min in exact arithmetic)
"""
import sys, os, itertools
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sympy as sp
from kernel_probe import positions_along_word, triangulation, detect_flip
from operators import classify
from float_kernel import gen_words


def flip_sequence(n, word, steps=600):
    """Discrete (i,j,k,l,removed) flip sequence of a braid word (geometry only locates the flips)."""
    labels = list(range(1, n + 4))
    cs = positions_along_word(n, word, steps)
    init = triangulation(cs[0], labels)
    prev, seq = init, []
    for fr in cs[1:]:
        cur = triangulation(fr, labels)
        if cur == prev:
            continue
        fl = detect_flip(prev, cur)
        if fl is None:
            prev = cur
            continue
        quad, removed = fl
        i, j, k, l, _ = classify(quad, {x: float(x) for x in labels}, fr)
        seq.append((i, j, k, l, frozenset(removed)))
        prev = cur
    return init, prev, seq


def _rblock(i, j, k, l):
    z = lambda a, b: sp.Integer(a) - sp.Integer(b)
    return sp.Matrix([[z(i, l) / z(i, k), z(i, j) / z(i, k)],
                      [z(l, k) / z(i, k), z(j, k) / z(i, k)]])


def _apply(state, i, j, k, l, removed):
    ijk, ikl = frozenset({i, j, k}), frozenset({i, k, l})
    ijl, jkl = frozenset({i, j, l}), frozenset({j, k, l})
    M = _rblock(i, j, k, l)
    out = dict(state)
    if set(removed) == {i, k}:
        a, b = out.pop(ijk, sp.Integer(0)), out.pop(ikl, sp.Integer(0))
        out[ijl] = out.get(ijl, 0) + M[0, 0] * a + M[0, 1] * b
        out[jkl] = out.get(jkl, 0) + M[1, 0] * a + M[1, 1] * b
    else:
        Mi = M.inv()
        a, b = out.pop(ijl, sp.Integer(0)), out.pop(jkl, sp.Integer(0))
        out[ijk] = out.get(ijk, 0) + Mi[0, 0] * a + Mi[0, 1] * b
        out[ikl] = out.get(ikl, 0) + Mi[1, 0] * a + Mi[1, 1] * b
    return {t: v for t, v in out.items() if v != 0}


def exact_generator(n, word):
    init, final, seq = flip_sequence(n, word)
    assert final == init, "braid did not return to initial triangulation"
    basis = sorted(init, key=lambda t: tuple(sorted(t)))
    idx = {t: r for r, t in enumerate(basis)}
    cols = []
    for t0 in basis:
        st = {t0: sp.Integer(1)}
        for (i, j, k, l, rem) in seq:
            st = _apply(st, i, j, k, l, rem)
        col = sp.zeros(len(basis), 1)
        for t, co in st.items():
            col[idx[t]] = sp.nsimplify(co)
        cols.append(col)
    return sp.Matrix.hstack(*cols)


def verify(n):
    I = sp.eye(2 * n + 1)
    N = {ij: exact_generator(n, w) - I for ij, w in gen_words(n).items()}
    ranks = all((Nm * Nm) == sp.zeros(*Nm.shape) and Nm.rank() == 2 for Nm in N.values())
    prods = all((N[a] * N[b]) == sp.zeros(*I.shape) for a, b in itertools.combinations(N, 2))
    comm = all((N[a] * N[b] - N[b] * N[a]) == sp.zeros(*I.shape) for a, b in itertools.combinations(N, 2))
    Mst = sp.Matrix([[N[ij][r, c] for r in range(2 * n + 1) for c in range(2 * n + 1)] for ij in N])
    indep, cn2 = Mst.rank(), n * (n - 1) // 2
    print(f"[n={n}] EXACT: rank-2 & N^2=0: {ranks}; all products N_ij N_kl=0: {prods}; "
          f"commute: {comm}; independent {indep}/{cn2}=C(n,2)")
    assert ranks, f"n={n}: a generator is not an exact rank-2 N^2=0 unipotent"
    assert prods and comm, f"n={n}: image is not exactly abelian"
    assert indep == cn2, f"n={n}: generators not independent ({indep} != {cn2})"
    return True


def main():
    print("== EXACT abelianization of inner PB_n (sympy; no floating point) ==\n")
    verify(3)
    verify(4)
    print("\n=> conditional-exact for inner PB_3, PB_4 (given the Delaunay-located flip sequences):")
    print("   rational f = abelianization, ker f = [PB_n, PB_n]; ker f' OPEN (conditional on C8' branch).")
    print("   General n and the full PB_{n+3} remain open (conjectural).")


if __name__ == "__main__":
    main()

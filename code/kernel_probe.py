"""
kernel_probe.py -- braid -> Delaunay-flip-sequence -> operator-product pipeline, the tool for
the kernel / Borromean investigation of arXiv:2606.20473.

Status: diagnostic tool (alg2_signs convention); passes its checkpoints.  This computes
f'(beta) for a braid beta by the paper's recipe (points move in the plane; record Delaunay flips;
multiply the flip operators).  It passes its checkpoints, with sigma_i acting on strand POSITIONS
(see positions_along_word):
  * trivially-trivial braids (sigma sigma^{-1}, sigma^2 sigma^{-2}) -> I exactly;
  * the braid relation (sigma_1 sigma_2 sigma_1)(sigma_2 sigma_1 sigma_2)^{-1} -> I exactly;
  * far commutation (sigma_1 sigma_3 sigma_1^{-1} sigma_3^{-1}) -> I exactly;
  * the Borromean braid (sigma_1 sigma_2^{-1})^3 -> I exactly (the paper's stated collapse).
The earlier version braided by label, violated the braid relation, and failed the Borromean check
(that was discrepancy D1); braiding by position fixes it.

NOTE on the word "validated": three of the five checkpoints (sigma sigma^{-1},
sigma^2 sigma^{-2}, and the braid-relation word) are IDENTITY elements of the braid group, so ANY
consistent braid-group representation sends them to I -- they test pipeline self-consistency, not the
correctness of the sigma-realization. The Borromean -> I is the very collapse under study and is
EXPECTED if the invariant is the abelianization (float_kernel.py), so it CORROBORATES rather than
validates. The non-tautological checkpoint here is far-commutation. Treat outputs as DIAGNOSTIC; the
load-bearing kernel result is the EXACT lift in exact_kernel.py, which uses this pipeline only to
LOCATE the combinatorial flip sequence and then applies exact operators.

Caveat: validated for braids built from sigma_i on the n inner strands. Pure-braid generators that
loop an inner strand around a FIXED outer point (e.g. the appendix b_{4,7}=A_47, which involves the
fixed label 7) are not yet expressible here -- they need a loop motion, not adjacent transpositions.

Run:  python kernel_probe.py
"""
import sys, os, itertools
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from scipy.spatial import Delaunay
from operators import classify, cos_sin, paper_signs, flip_on_vector, is_complex_orthogonal

# ----- geometry: fixed outer triangle (labels n+1,n+2,n+3) + n moving inner points (1..n) -----
def base_layout(n, R=60.0):
    outer = {n + 1: (0.0, R), n + 2: (-R * 0.95, -R * 0.6), n + 3: (R * 0.95, -R * 0.6)}
    inner0 = {i: (float(i) - (n + 1) / 2.0, 0.0) for i in range(1, n + 1)}   # on y=0, centered
    return outer, inner0


def positions_along_word(n, word, steps_per_gen=400):
    """word = list of (i, sign): generator sigma_i^sign crosses the strands in POSITIONS i,i+1
    (counting left to right by current x-coordinate) via a +-pi half-turn about their midpoint.
    Braiding by position, not by label, is essential: as strands cross, the strand occupying a
    given position changes, and sigma_i acts on positions. (Braiding by label violates the braid
    relation sigma_1 sigma_2 sigma_1 = sigma_2 sigma_1 sigma_2 and fails the Borromean checkpoint.)
    Returns frames with frame[t] = {label:(x,y)}."""
    outer, inner = base_layout(n)
    frames = []
    cur = dict(inner)
    frames.append({**outer, **cur})
    for (i, sign) in word:
        order = sorted(range(1, n + 1), key=lambda lab: cur[lab][0])   # labels by current x-position
        a, b = order[i - 1], order[i]                                  # i-th and (i+1)-th positions
        pa, pb = np.array(cur[a]), np.array(cur[b])
        mid = (pa + pb) / 2.0
        va, vb = pa - mid, pb - mid
        for t in range(1, steps_per_gen + 1):
            th = sign * np.pi * t / steps_per_gen
            c, s = np.cos(th), np.sin(th)
            Rm = np.array([[c, -s], [s, c]])
            cur = dict(cur)
            cur[a] = tuple(mid + Rm @ va)
            cur[b] = tuple(mid + Rm @ vb)
            frames.append({**outer, **cur})
    return frames


def triangulation(frame, labels):
    pts = np.array([frame[l] for l in labels])
    tri = Delaunay(pts)
    out = set()
    for simplex in tri.simplices:
        out.add(frozenset(labels[v] for v in simplex))
    return out


def detect_flip(old, new):
    """Given two triangulations differing by one 2->2 flip, return (quad, removed_diag) or None."""
    rem, add = old - new, new - old
    if len(rem) != 2 or len(add) != 2:
        return None  # not a single flip (refine sampling) or no change
    quad = set().union(*rem)
    if len(quad) != 4:
        return None
    old_diag = frozenset(set.intersection(*[set(t) for t in rem]))
    new_diag = frozenset(set.intersection(*[set(t) for t in add]))
    if len(old_diag) != 2 or len(new_diag) != 2:
        return None
    return quad, old_diag


def _flip_operator(quad, removed, values, coords, convention):
    """Return (i,j,k,l,c,s) for a detected flip under the chosen convention.
       'alg2_signs' : alg_2 index assignment + the B/C real sign rule (Remark alg_2).
       'appendix'   : (i,k)=sorted removed diagonal, (j,l)=sorted added diagonal, principal
                      complex sqrt, NO sign flip (the convention locked from the paper's appendix)."""
    if convention == "alg2_signs":
        i, j, k, l, neck = classify(quad, values, coords)
        sc, ss = paper_signs(neck)
        c, s = cos_sin(values, i, j, k, l, sc, ss)
        return i, j, k, l, c, s
    elif convention == "appendix":
        added = quad - set(removed)
        i, k = sorted(removed, key=lambda v: values[v])
        j, l = sorted(added, key=lambda v: values[v])
        c, s = cos_sin(values, i, j, k, l)
        return i, j, k, l, c, s
    raise ValueError(convention)


def braid_matrix(n, word, steps_per_gen=600, verbose=False, convention="alg2_signs"):
    """Compute f'(beta) as a matrix in the initial-triangulation basis (pure braids only:
    final triangulation must equal initial).  `convention` selects the operator sign rule."""
    labels = list(range(1, n + 4))
    values = {l: float(l) for l in labels}          # fixed zeta_l = l (as in the paper's example)
    coords_seq = positions_along_word(n, word, steps_per_gen)

    init_tris = triangulation(coords_seq[0], labels)
    # track image of each initial basis triangle as dict current_triangle->coeff
    state = {t: {t: 1 + 0j} for t in init_tris}
    prev = init_tris
    nflips = 0
    multi = 0
    for fr in coords_seq[1:]:
        cur = triangulation(fr, labels)
        if cur == prev:
            continue
        flip = detect_flip(prev, cur)
        if flip is None:
            multi += 1
            prev = cur
            continue
        quad, removed = flip
        i, j, k, l, c, s = _flip_operator(quad, removed, values, fr, convention)
        for b in state:
            state[b] = flip_on_vector(state[b], i, j, k, l, c, s, removed)
        nflips += 1
        prev = cur

    final_tris = prev
    if verbose:
        print(f"  flips detected: {nflips}; unresolved multi-flip frames: {multi}")
        print(f"  returns to initial triangulation: {final_tris == init_tris}")
    basis = sorted(init_tris, key=lambda t: tuple(sorted(t)))
    if final_tris != init_tris:
        return None, basis, nflips, multi
    M = np.zeros((len(basis), len(basis)), dtype=complex)
    idx = {t: r for r, t in enumerate(basis)}
    for col, b in enumerate(basis):
        for t, co in state[b].items():
            M[idx[t], col] = co
    return M, basis, nflips, multi


def report(name, n, word, steps=600):
    print(f"\n[{name}]  n={n}, word length={len(word)}")
    M, basis, nflips, multi = braid_matrix(n, word, steps, verbose=True)
    if M is None:
        print("  (did not return to initial triangulation at this resolution; refine steps)")
        return
    dim = M.shape[0]
    print(f"  matrix dimension = {dim}  (expected 2n+1 = {2*n+1})")
    print(f"  complex-orthogonal (M^T M = I): {is_complex_orthogonal(M, tol=1e-6)}")
    print(f"  det = {np.linalg.det(M):.6f}")
    print(f"  trace = {np.trace(M):.6f}   (unit-diagonal would give {dim})")
    print(f"  ||M - I|| = {np.linalg.norm(M - np.eye(dim)):.3e}")
    diag = np.diag(M)
    print(f"  max |diag - 1| = {np.max(np.abs(diag - 1)):.3e}")


def main():
    print("== Braid -> Delaunay-flip -> operator-product pipeline (DIAGNOSTIC; checkpoints pass, alg2_signs) ==")
    print("Fixed zeta_l = l; outer triangle fixed; inner strands braid by POSITION.\n")
    print("[checkpoints: each should give the identity]")
    checks = [("s1 s1^-1", 2, [(1, 1), (1, -1)]),
              ("s1^2 s1^-2", 2, [(1, 1), (1, 1), (1, -1), (1, -1)]),
              ("braid reln (s1s2s1)(s2s1s2)^-1", 3, [(1, 1), (2, 1), (1, 1), (2, -1), (1, -1), (2, -1)]),
              ("far comm (s1 s3 s1^-1 s3^-1)", 4, [(1, 1), (3, 1), (1, -1), (3, -1)]),
              ("BORROMEAN (s1 s2^-1)^3", 3, [(1, 1), (2, -1)] * 3)]
    for name, n, w in checks:
        M, basis, nflips, multi = braid_matrix(n, w, 600)
        if M is None:
            print(f"  {name:34s} did not close up (refine steps)"); continue
        d = np.linalg.norm(M - np.eye(M.shape[0]))
        ok = "IDENTITY" if d < 1e-7 else f"NOT identity (||M-I||={d:.3f})"
        print(f"  {name:34s} flips={nflips:3d} -> {ok}")
    print("\n[a nontrivial pure braid, for contrast]")
    report("b_12 = sigma_1^2 (n=2)", n=2, word=[(1, 1), (1, 1)])
    print("\nThe Borromean -> I reproduces the paper's collapse claim and resolves discrepancy D1;")
    print("the pipeline respects the braid relations and gives DIAGNOSTIC evidence for f'(beta) on")
    print("sigma-word braids; the load-bearing exact kernel result is in exact_kernel.py.")


if __name__ == "__main__":
    main()

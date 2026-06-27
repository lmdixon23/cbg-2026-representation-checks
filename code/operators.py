"""
operators.py -- the complex-orthogonal flip operators of the colored-braid-groupoid
representation of Rohozhkin, "Invariants of the Colored Braid Groupoid" (arXiv:2606.20473).

Everything here is built from the paper's formulas (eq. operator-2/3 and the quadrilateral
sign rule of Remark "alg_2" / Fig. "quadrilaterials"), and is exercised by relations.py.

Conventions
-----------
* A "configuration" is a dict  label -> real value zeta_label  (all distinct), plus a dict
  label -> (x, y) planar coordinate.  zeta plays the role of the paper's variable; the planar
  coordinates fix the Delaunay/flip geometry.
* A quadrilateral is a 4-set of labels in convex position.  A flip swaps one diagonal for the
  other.  Triangles are identified by frozenset-of-3-labels (vertex order is irrelevant).
* cos/sin use the paper's square-root formulas; over C they are taken on the principal branch,
  and the quadrilateral-type sign rule (below) selects the correct branch.

Key structural fact (Lemma in the accompanying note): cos^2 + sin^2 = 1 is *equivalent* to the
Ptolemy/Plucker identity  z_il z_jk + z_ij z_kl = z_ik z_jl , which holds identically, so every
flip operator is complex-orthogonal (A^T A = I) for ANY ordering of the points -- it is
genuine orthogonality only when the ratios are positive, complex-orthogonality otherwise.
"""
from __future__ import annotations
import cmath
import numpy as np


# ----------------------------------------------------------------------------------
# trig functions and the Ptolemy identity
# ----------------------------------------------------------------------------------
def zdiff(values, a, b):
    return values[a] - values[b]


def cos_sin(values, i, j, k, l, sign_c=1, sign_s=1):
    """Paper's cos/sin for the flip ik->jl, with optional branch signs (+-1).
    Returns complex (cos phi, sin phi); principal branch of sqrt."""
    z = lambda a, b: values[a] - values[b]
    c = sign_c * cmath.sqrt((z(i, l) * z(j, k)) / (z(i, k) * z(j, l)))
    s = sign_s * cmath.sqrt((z(i, j) * z(k, l)) / (z(i, k) * z(j, l)))
    return c, s


def ptolemy_residual(values, i, j, k, l):
    """z_il z_jk + z_ij z_kl - z_ik z_jl.  Identically 0  ==>  cos^2 + sin^2 = 1."""
    z = lambda a, b: values[a] - values[b]
    return z(i, l) * z(j, k) + z(i, j) * z(k, l) - z(i, k) * z(j, l)


# ----------------------------------------------------------------------------------
# geometry: convex CCW order and quadrilateral classification
# ----------------------------------------------------------------------------------
def ccw_cycle(quad, coords):
    """The 4 labels of `quad` in convex counter-clockwise order (by angle about centroid)."""
    cx = sum(coords[v][0] for v in quad) / 4.0
    cy = sum(coords[v][1] for v in quad) / 4.0
    return sorted(quad, key=lambda v: np.arctan2(coords[v][1] - cy, coords[v][0] - cx))


def classify(quad, values, coords):
    """Implements Remark alg_2.  Returns (i, j, k, l, necklace) where:
       i = min-zeta vertex; k = its cyclic-opposite (so {i,k} is a diagonal);
       j,l = the other two with zeta_j < zeta_l;
       necklace = ranks (1..4 by zeta) read CCW starting at i  -- the quadrilateral "type"."""
    cyc = ccw_cycle(quad, coords)
    i = min(quad, key=lambda v: values[v])
    p = cyc.index(i)
    k = cyc[(p + 2) % 4]
    rest = [v for v in quad if v not in (i, k)]
    j, l = sorted(rest, key=lambda v: values[v])
    order = sorted(quad, key=lambda v: values[v])
    rank = {v: r + 1 for r, v in enumerate(order)}
    necklace = tuple(rank[v] for v in (cyc[p:] + cyc[:p]))
    return i, j, k, l, necklace


# The paper's sign rule (Remark alg_2 / Fig. quadrilaterials): s=1 (cos flip) for type B,
# r=1 (sin flip) for type C, else 0.  In the CCW-necklace encoding used here:
#    type B = (1,3,2,4)   ->  cos sign -1
#    type C = (1,2,4,3)   ->  sin sign -1
# These two assignments were singled out by an exhaustive search over all 4^6 type->sign
# rules (sign_uniqueness.py): exactly two rules make the pentagon hold for all 120 orderings;
# this one IS the paper's, and the other is its orientation mirror.
_COS_MINUS_TYPE = (1, 3, 2, 4)   # "B"
_SIN_MINUS_TYPE = (1, 2, 4, 3)   # "C"


def paper_signs(necklace):
    """(sign_c, sign_s) in {+1,-1}^2 per the paper's quadrilateral-type rule."""
    return (-1 if necklace == _COS_MINUS_TYPE else 1,
            -1 if necklace == _SIN_MINUS_TYPE else 1)


# ----------------------------------------------------------------------------------
# flip operators
# ----------------------------------------------------------------------------------
def flip_on_vector(vec, i, j, k, l, c, s, removed_diagonal):
    """Apply the flip operator to a vector (dict: triangle frozenset -> complex coeff).

    The operator (paper eq. operator-2/3) sends, for the diagonal {i,k} -> {j,l}:
        D_ijk |-> c D_ijl - s D_jkl ,   D_ikl |-> s D_ijl + c D_jkl ,
    and fixes every other triangle.  If the actual flip removes {j,l} instead (i.e. it is the
    inverse flip jl->ik), the inverse operator is applied.  Triangles are frozensets."""
    ijk, ikl = frozenset({i, j, k}), frozenset({i, k, l})
    ijl, jkl = frozenset({i, j, l}), frozenset({j, k, l})
    forward = (set(removed_diagonal) == {i, k})
    out = {}
    for t, co in vec.items():
        if forward and t == ijk:
            out[ijl] = out.get(ijl, 0) + co * c
            out[jkl] = out.get(jkl, 0) - co * s
        elif forward and t == ikl:
            out[ijl] = out.get(ijl, 0) + co * s
            out[jkl] = out.get(jkl, 0) + co * c
        elif (not forward) and t == ijl:                       # inverse: {j,l} -> {i,k}
            out[ijk] = out.get(ijk, 0) + co * c
            out[ikl] = out.get(ikl, 0) + co * s
        elif (not forward) and t == jkl:
            out[ijk] = out.get(ijk, 0) - co * s
            out[ikl] = out.get(ikl, 0) + co * c
        else:
            out[t] = out.get(t, 0) + co
    return out


def givens_2x2(c, s):
    """The active 2x2 block [[c, s], [-s, c]] (complex-orthogonal: block^T block = I when c^2+s^2=1)."""
    return np.array([[c, s], [-s, c]], dtype=complex)


def is_complex_orthogonal(M, tol=1e-9):
    """A^T A = I (NOT conjugate transpose) -- the relevant notion here."""
    M = np.asarray(M, dtype=complex)
    return np.linalg.norm(M.T @ M - np.eye(M.shape[0])) < tol


if __name__ == "__main__":
    # smoke test: Ptolemy identity and orthogonality of a single flip block
    import itertools
    vals = {1: 5.0, 2: 4.0, 3: 3.0, 4: 2.0}
    worst = max(abs(ptolemy_residual(vals, *p)) for p in itertools.permutations([1, 2, 3, 4]))
    c, s = cos_sin(vals, 1, 2, 3, 4)
    print(f"max |Ptolemy residual| over 4! orderings = {worst:.2e}")
    print(f"single block complex-orthogonal: {is_complex_orthogonal(givens_2x2(c, s))}")

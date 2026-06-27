"""
appendix_full_matrix.py -- full-matrix check of the paper's two worked pure braids, and an erratum.

The paper (arXiv:2606.20473) prints, for zeta_l = l, the complete 9x9 matrices f'(b_{4,7}) and
f'(b_{5,6}) BOTH as an explicit product of named complex-orthogonal flip matrices AND as a single
displayed "= (9x9)" result. This script reconstructs the printed factors and multiplies them, and
audits the displayed results. Two facts, both COMPUTED here (no claim is printed without a check):

(1) The paper's FACTORS are correct. Multiplying the printed flip matrices in printed (left-to-right)
    order yields, in both examples, a matrix that IS complex-orthogonal, has determinant 1, is
    unit-diagonal, and is a rank-2 one-step unipotent ((M-I)^2 = 0). This matches the paper's stated
    properties and the rational computation in appendix_reconstruction.py via f' = D f D^{-1}
    (note, Prop. "diagonally similar"). The rank-2 unipotent form independently corroborates the
    abelianization structure (float_kernel.py) on the inner-outer generators A_47 and A_56.

(2) ERRATUM: the two DISPLAYED result matrices are wrong. Each is unit-diagonal but equals
    I + (skew-symmetric), hence is NOT complex-orthogonal and has determinant 617.52 and 49
    respectively -- contradicting the paper's own properties AND the product of its own factors.

Factors are stored as sparse overrides on the 9x9 identity: each entry (r, c, v) overwrites I[r,c].
'v' ending in 'j' is imaginary; v == '0' zeroes a (formerly unit) diagonal entry.

Run:  python appendix_full_matrix.py
"""
import os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

def _v(s):
    return complex(0.0, float(s[:-1])) if s.endswith("j") else complex(float(s), 0.0)

def _mat(overrides):
    M = np.eye(9, dtype=complex)
    for r, c, v in overrides:
        M[r, c] = _v(v)
    return M

# --- f'(b_{4,7}) = A_4615 A_1524 A_4716 A_3417 A_1645 A_1746 A_2437 A_3714 A_4627 A_2734 A_1467 A_6724 A_4516 A_1647 A_2415 A_1546
B47_FACTORS = [
    [(2,2,'0.790569'),(2,8,'-0.612372'),(3,2,'0.612372'),(3,3,'0'),(3,8,'0.790569'),(4,3,'1'),(4,4,'0'),(5,4,'1'),(5,5,'0'),(6,5,'1'),(6,6,'0'),(7,6,'1'),(7,7,'0'),(8,7,'1'),(8,8,'0')],
    [(0,0,'1.06066'),(0,2,'-0.353553j'),(2,2,'0'),(2,3,'1'),(3,3,'0'),(3,4,'1'),(4,4,'0'),(4,5,'1'),(5,0,'0.353553j'),(5,2,'1.06066'),(5,5,'0')],
    [(3,3,'0.894427'),(3,8,'-0.447214'),(4,3,'0.447214'),(4,4,'0'),(4,8,'0.894427'),(5,4,'1'),(5,5,'0'),(6,5,'1'),(6,6,'0'),(7,6,'1'),(7,7,'0'),(8,7,'1'),(8,8,'0')],
    [(1,1,'1.41421'),(1,6,'-1j'),(3,1,'1j'),(3,3,'0'),(3,6,'1.41421'),(4,3,'1'),(4,4,'0'),(5,4,'1'),(5,5,'0'),(6,5,'1'),(6,6,'0')],
    [(2,2,'1.26491'),(2,3,'0.774597j'),(3,3,'0'),(3,4,'1'),(4,4,'0'),(4,5,'1'),(5,5,'0'),(5,6,'1'),(6,6,'0'),(6,7,'1'),(7,2,'-0.774597j'),(7,3,'1.26491'),(7,7,'0')],
    [(2,2,'1.11803'),(2,4,'0.5j'),(4,4,'0'),(4,5,'1'),(5,5,'0'),(5,6,'1'),(6,6,'0'),(6,7,'1'),(7,7,'0'),(7,8,'1'),(8,2,'-0.5j'),(8,4,'1.11803'),(8,8,'0')],
    [(5,5,'0.790569'),(5,6,'0.612372'),(6,6,'0'),(6,7,'1'),(7,7,'0'),(7,8,'1'),(8,5,'-0.612372'),(8,6,'0.790569'),(8,8,'0')],
    [(1,1,'0.707107'),(1,8,'-0.707107'),(2,1,'0.707107'),(2,2,'0'),(2,8,'0.707107'),(3,2,'1'),(3,3,'0'),(4,3,'1'),(4,4,'0'),(5,4,'1'),(5,5,'0'),(6,5,'1'),(6,6,'0'),(7,6,'1'),(7,7,'0'),(8,7,'1'),(8,8,'0')],
    [(5,5,'1.09545'),(5,8,'0.447214j'),(7,5,'-0.447214j'),(7,7,'0'),(7,8,'1.09545'),(8,7,'1'),(8,8,'0')],
    [(4,4,'1.26491'),(4,6,'-0.774597j'),(6,6,'0'),(6,7,'1'),(7,4,'0.774597j'),(7,6,'1.26491'),(7,7,'0')],
    [(2,2,'0'),(2,4,'1'),(3,2,'-2j'),(3,3,'2.23607'),(4,4,'0'),(4,5,'1'),(5,5,'0'),(5,6,'1'),(6,6,'0'),(6,7,'1'),(7,7,'0'),(7,8,'1'),(8,2,'-2.23607'),(8,3,'-2j'),(8,8,'0')],
    [(6,6,'0'),(6,7,'2.23607j'),(6,8,'-2.44949'),(7,7,'2.44949'),(7,8,'2.23607j'),(8,6,'1'),(8,8,'0')],
    [(2,2,'1.26491'),(2,7,'-0.774597j'),(4,2,'0.774597j'),(4,4,'0'),(4,7,'1.26491'),(5,4,'1'),(5,5,'0'),(6,5,'1'),(6,6,'0'),(7,6,'1'),(7,7,'0')],
    [(3,3,'0.894427'),(3,4,'0.447214'),(4,4,'0'),(4,5,'1'),(5,5,'0'),(5,6,'1'),(6,6,'0'),(6,7,'1'),(7,7,'0'),(7,8,'1'),(8,3,'-0.447214'),(8,4,'0.894427'),(8,8,'0')],
    [(0,0,'1.06066'),(0,5,'0.353553j'),(2,0,'-0.353553j'),(2,2,'0'),(2,5,'1.06066'),(3,2,'1'),(3,3,'0'),(4,3,'1'),(4,4,'0'),(5,4,'1'),(5,5,'0')],
    [(2,2,'0.790569'),(2,3,'0.612372'),(3,3,'0'),(3,4,'1'),(4,4,'0'),(4,5,'1'),(5,5,'0'),(5,6,'1'),(6,6,'0'),(6,7,'1'),(7,7,'0'),(7,8,'1'),(8,2,'-0.612372'),(8,3,'0.790569'),(8,8,'0')],
]
# --- f'(b_{5,6}) = A_4625 A_5716 A_1645 A_2567 A_6715 A_4526 A_1546 A_2657
B56_FACTORS = [
    [(6,6,'0.816497'),(6,8,'-0.57735'),(7,6,'0.57735'),(7,7,'0'),(7,8,'0.816497'),(8,7,'1'),(8,8,'0')],
    [(3,3,'0.774597'),(3,8,'-0.632456'),(4,3,'0.632456'),(4,4,'0'),(4,8,'0.774597'),(5,4,'1'),(5,5,'0'),(6,5,'1'),(6,6,'0'),(7,6,'1'),(7,7,'0'),(8,7,'1'),(8,8,'0')],
    [(2,2,'1.26491'),(2,3,'-0.774597j'),(3,3,'0'),(3,4,'1'),(4,4,'0'),(4,5,'1'),(5,5,'0'),(5,6,'1'),(6,6,'0'),(6,7,'1'),(7,2,'0.774597j'),(7,3,'1.26491'),(7,7,'0')],
    [(7,7,'1.29099j'),(7,8,'1.63299'),(8,7,'-1.63299'),(8,8,'1.29099j')],
    [(3,3,'-1.22474j'),(3,8,'-1.58114'),(4,3,'1.58114'),(4,4,'0'),(4,8,'-1.22474j'),(5,4,'1'),(5,5,'0'),(6,5,'1'),(6,6,'0'),(7,6,'1'),(7,7,'0'),(8,7,'1'),(8,8,'0')],
    [(5,5,'1.22474'),(5,7,'-0.707107j'),(6,5,'0.707107j'),(6,6,'0'),(6,7,'1.22474'),(7,6,'1'),(7,7,'0')],
    [(2,2,'0.790569'),(2,3,'0.612372'),(3,3,'0'),(3,4,'1'),(4,4,'0'),(4,5,'1'),(5,5,'0'),(5,6,'1'),(6,6,'0'),(6,7,'1'),(7,2,'-0.612372'),(7,3,'0.790569'),(7,7,'0')],
    [(7,7,'0.790569'),(7,8,'0.612372'),(8,7,'-0.612372'),(8,8,'0.790569')],
]

B47_PRINTED = np.array([
    [1,-0.35355,0,0,1.11803,0.27386,0,0,-1.09545],[0.35355,1,-1,0,0,0,-1.06066,0,0],
    [0,1,1,0,-3.16228,-0.7746,0,0,3.09839],[0,0,0,1,0,0,0,0,0],
    [-1.11803,0,3.16228,0,1,0,-3.3541,0,0],[-0.27386,0,0.7746,0,0,1,-0.82158,0,0],
    [0,1.06066,0,0,3.3541,0.82158,1,0,3.28634],[0,0,0,0,0,0,0,1,0],
    [1.09545,0,-3.09839,0,0,0,-3.28634,0,1]], dtype=complex)
B56_PRINTED = np.array([
    [1,0,0,0,0,0,0,0,0],[0,1,0,0,0,0,0,0,0],[0,0,1,0.7746,-0.94868,0,0,-0.75,0.96825],
    [0,0,-0.7746,1,1.22474,0,0.7303,-0.06455,-1.25],[0,0,0.94868,-1.22474,1,0,-0.89443,1.26491,0],
    [0,0,0,0,0,1,0,0,0],[0,0,0,-0.7303,0.89443,0,1,-0.70711,0.91287],
    [0,0,0.75,0.06455,-1.26491,0,0.70711,1,-1.29099],[0,0,-0.96825,1.25,0,0,-0.91287,1.29099,1]], dtype=complex)


def product(factors):
    M = np.eye(9, dtype=complex)
    for ov in factors:
        M = M @ _mat(ov)          # printed (left-to-right) order
    return M

def structure(M):
    I = np.eye(9); N = M - I
    return dict(unit=np.allclose(np.diag(M), 1, atol=1e-4),
               det=np.linalg.det(M),
               orth=np.allclose(M.T @ M, I, atol=1e-3),
               rank=np.linalg.matrix_rank(N, tol=1e-4),
               nilp=np.allclose(N @ N, 0, atol=1e-3),
               skew=np.allclose(N + N.T, 0, atol=1e-3))


def rational_check():
    """Transcription-INDEPENDENT witness of the correct value. The rational f(b), built straight from
    the operator formulas (no printed digits), has the same det/rank/unipotency as the correct f' by
    the diagonal similarity f'=D f D^{-1} (C8). So the correct f'(b) is det-1, rank-2 unipotent --
    which the displayed result matrices (det 617.52 / 49, rank 4) are not, independent of any factor
    transcription."""
    from appendix_reconstruction import full_product
    out = {}
    for name in ("rational",):
        M = full_product("rational", "reversed")   # rational f(b_4,7) from operator formulas only
        N = M - np.eye(M.shape[0])
        out = dict(det=np.linalg.det(M).real,
                   rank=np.linalg.matrix_rank(N, tol=1e-9),
                   nilp=np.allclose(N @ N, 0, atol=1e-9))
    return out


def main():
    print("== Appendix worked examples: full-matrix audit (arXiv:2606.20473) ==\n")
    ok = True
    for name, facs, printed in [("b_4,7", B47_FACTORS, B47_PRINTED), ("b_5,6", B56_FACTORS, B56_PRINTED)]:
        P = product(facs); sp = structure(P); sd = structure(printed)
        print(f"f'({name}):")
        print(f"  product of {len(facs)} printed factors  : orthogonal={sp['orth']}, det={sp['det'].real:.4f}, "
              f"unit-diag={sp['unit']}, rank(M-I)={sp['rank']}, (M-I)^2=0={sp['nilp']}  <- correct, rank-2 unipotent")
        print(f"  displayed result matrix (as printed): orthogonal={sd['orth']}, det={sd['det'].real:.4f}, "
              f"unit-diag={sd['unit']}, rank(M-I)={sd['rank']}, I+skew={sd['skew']}  <- ERRATUM")
        differ = not np.allclose(P, printed, atol=1e-2)
        print(f"  product != displayed result: {differ}\n")
        # hard gates
        assert sp['orth'] and abs(sp['det'] - 1) < 1e-2 and sp['unit'] and sp['rank'] == 2 and sp['nilp'], \
            f"{name}: factor product is not the expected orthogonal det-1 rank-2 unipotent"
        assert (not sd['orth']) and abs(sd['det'] - 1) > 1e-2, f"{name}: displayed matrix unexpectedly consistent"
        assert differ, f"{name}: displayed matrix unexpectedly equals the factor product"
        ok = ok and differ
    # transcription-independent corroboration of the correct value (b_4,7)
    rc = rational_check()
    print(f"transcription-INDEPENDENT check (rational f(b_4,7) from formulas, no printed digits): "
          f"det={rc['det']:.4f}, rank(M-I)={rc['rank']}, (M-I)^2=0={rc['nilp']}")
    assert abs(rc['det'] - 1) < 1e-9 and rc['rank'] == 2 and rc['nilp'], "rational witness failed"
    print("  => the correct f'(b) is det-1 rank-2 unipotent regardless of the printed factor digits;")
    print("     the displayed matrices (det 617.52 / 49, rank 4) cannot be it. Erratum stands independently.\n")
    print("RESULT: in both examples the product of the paper's own factors is complex-orthogonal,")
    print("det 1, unit-diagonal, rank-2 unipotent; the displayed result matrices are erroneous (I+skew).")
    print("The rank-2 unipotent factor products corroborate the abelianization on A_47, A_56.")


if __name__ == "__main__":
    main()

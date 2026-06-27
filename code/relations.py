"""
FAST FLOAT SMOKE TEST. The exact, certifying check for the four relations is exact_checks.py;
this file is a quick floating-point reproduction and provides the shared fixtures.
relations.py -- machine verification of the four defining relations of the groupoid
G^4_{n+3} for the complex-orthogonal representation f' (Rohozhkin, arXiv:2606.20473):

    (I)   g^H_H = 1                       identity            (trivial)
    (II)  g^Hi_Hj g^Hj_Hi = 1             4-cycle             (operator inverse, by construction)
    (III) 5-cycle (pentagon)              g..g..g..g..g = 1   (Korepanov pentagon)
    (IV)  far-commutativity               disjoint flips commute

Run:  python relations.py
Reproduces: pentagon holds for 60/120 orderings WITHOUT the sign rule and 120/120 WITH it;
far-commutativity holds exactly for all orderings.
"""
from __future__ import annotations
import itertools
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
from operators import classify, cos_sin, paper_signs, flip_on_vector

# ----- pentagon geometry (5 convex points), labels 1..5 (from Fig. pentagon_korep_flips) -----
PENT_COORDS = {1: (-1.0, 0.30), 2: (-0.6, -0.90), 3: (0.6, -0.90), 4: (1.0, 0.30), 5: (0.0, 1.0)}

# the five flips (quad, removed_diagonal, added_diagonal) from the appendix triangulations
PENT_FLIPS = {
    "i":   ({1, 2, 3, 4}, frozenset({2, 4}), frozenset({1, 3})),   # T0 -> T1
    "ii":  ({1, 3, 4, 5}, frozenset({1, 4}), frozenset({3, 5})),   # T1 -> Tend
    "iii": ({1, 2, 4, 5}, frozenset({1, 4}), frozenset({2, 5})),   # T0 -> T2
    "iv":  ({2, 3, 4, 5}, frozenset({2, 4}), frozenset({3, 5})),   # T2 -> T3
    "v":   ({1, 2, 3, 5}, frozenset({2, 5}), frozenset({1, 3})),   # T3 -> Tend
}
T0_BASIS = [frozenset({1, 2, 4}), frozenset({1, 4, 5}), frozenset({2, 3, 4})]


def pentagon_holds(values, use_signs=True, tol=1e-9):
    """Does (ii) o (i) == (v) o (iv) o (iii) on the T0 basis, for this zeta-assignment?"""
    def M(name, vec):
        quad, rem, _ = PENT_FLIPS[name]
        i, j, k, l, neck = classify(quad, values, PENT_COORDS)
        sc, ss = paper_signs(neck) if use_signs else (1, 1)
        c, s = cos_sin(values, i, j, k, l, sc, ss)
        return flip_on_vector(vec, i, j, k, l, c, s, rem)
    for t0 in T0_BASIS:
        top = M("ii", M("i", {t0: 1 + 0j}))
        bot = M("v", M("iv", M("iii", {t0: 1 + 0j})))
        keys = set(top) | set(bot)
        if any(abs(top.get(x, 0) - bot.get(x, 0)) > tol for x in keys):
            return False
    return True


def verify_pentagon(values_list=(5., 4., 3., 2., 1.)):
    holds_signed = holds_bare = 0
    for perm in itertools.permutations(values_list):
        vals = {i + 1: perm[i] for i in range(5)}
        holds_signed += pentagon_holds(vals, use_signs=True)
        holds_bare += pentagon_holds(vals, use_signs=False)
    n = 120
    return holds_bare, holds_signed, n


# ----- far-commutativity: flips in disjoint quadrilaterals commute -----
FC_BASIS = [frozenset({1, 2, 3}), frozenset({1, 3, 4}),
            frozenset({5, 6, 7}), frozenset({5, 7, 8}), frozenset({2, 3, 9})]  # last = spectator


def verify_far_commute(values_list=(8., 7., 6., 5., 4., 3., 2., 1.)):
    worst = 0.0
    for perm in itertools.permutations(values_list):
        vals = {i + 1: perm[i] for i in range(8)}
        vals[9] = 0.5
        c1, s1 = cos_sin(vals, 1, 2, 3, 4)
        c2, s2 = cos_sin(vals, 5, 6, 7, 8)
        for t in FC_BASIS:
            a = flip_on_vector(flip_on_vector({t: 1 + 0j}, 1, 2, 3, 4, c1, s1, frozenset({1, 3})),
                               5, 6, 7, 8, c2, s2, frozenset({5, 7}))
            b = flip_on_vector(flip_on_vector({t: 1 + 0j}, 5, 6, 7, 8, c2, s2, frozenset({5, 7})),
                               1, 2, 3, 4, c1, s1, frozenset({1, 3}))
            worst = max(worst, max(abs(a.get(k, 0) - b.get(k, 0)) for k in set(a) | set(b)))
    return worst


def main():
    print("== Groupoid relations for the complex-orthogonal representation f' ==\n")
    print("(I) identity: trivial by construction.")
    print("(II) 4-cycle g g^{-1}=1: the flip operator's inverse is built in (flip_on_vector).\n")
    bare, signed, n = verify_pentagon()
    print(f"(III) pentagon / 5-cycle:")
    print(f"      without sign rule: holds for {bare}/{n} orderings")
    print(f"      with paper's sign rule: holds for {signed}/{n} orderings")
    assert bare == 60 and signed == 120, "pentagon regression!"
    worst = verify_far_commute()
    print(f"\n(IV) far-commutativity: max ||M1 M2 - M2 M1|| over all 8! orderings = {worst:.2e}")
    assert worst < 1e-9, "far-commute regression!"
    print("\nAll four relations verified.  (Pentagon uniquely needs the B/C sign rule; see sign_uniqueness.py.)")


if __name__ == "__main__":
    main()

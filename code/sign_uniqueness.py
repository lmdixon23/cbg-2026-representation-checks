"""
FAST FLOAT SMOKE TEST. The exact, certifying check for sign-rule uniqueness is exact_checks.py;
this file is a quick floating-point reproduction and provides the shared fixtures.
sign_uniqueness.py -- exhaustive search establishing that the paper's quadrilateral-type sign
rule is essentially the unique fix that makes the pentagon hold for all point orderings.

Over all 4^6 = 4096 ways of attaching a (cos-sign, sin-sign) in {+-1}^2 to each of the six
quadrilateral types, exactly TWO make the 5-cycle relation hold for all 120 zeta-orderings.
Both have the form "cos-flip on one type, sin-flip on one (other) type"; one is the paper's
B/C rule, the other is its orientation mirror (necklaces are cyclic reverses).

Run:  python sign_uniqueness.py
"""
import sys, os, itertools
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from operators import classify, cos_sin, flip_on_vector
from relations import PENT_COORDS, PENT_FLIPS, T0_BASIS

CHOICES = list(itertools.product([1, -1], repeat=2))   # (cos-sign, sin-sign)


def collect():
    """For each ordering, precompute per-flip (i,j,k,l, base c0,s0, necklace, removed)."""
    types = {}
    data = []
    for perm in itertools.permutations([5., 4., 3., 2., 1.]):
        vals = {i + 1: perm[i] for i in range(5)}
        info = {}
        for name, (quad, rem, _) in PENT_FLIPS.items():
            i, j, k, l, neck = classify(quad, vals, PENT_COORDS)
            c0, s0 = cos_sin(vals, i, j, k, l)
            types.setdefault(neck, len(types))
            info[name] = (i, j, k, l, c0, s0, types[neck], rem)
        data.append(info)
    return data, types


def holds(info, sgn, tol=1e-9):
    def M(name, vec):
        i, j, k, l, c0, s0, ty, rem = info[name]
        sc, ss = sgn[ty]
        return flip_on_vector(vec, i, j, k, l, sc * c0, ss * s0, rem)
    for t0 in T0_BASIS:
        top = M("ii", M("i", {t0: 1 + 0j}))
        bot = M("v", M("iv", M("iii", {t0: 1 + 0j})))
        if any(abs(top.get(x, 0) - bot.get(x, 0)) > tol for x in set(top) | set(bot)):
            return False
    return True


def main():
    data, types = collect()
    nt = len(types)
    inv = {v: k for k, v in types.items()}
    print(f"{nt} quadrilateral types encountered (CCW rank-necklaces): {list(types)}")
    winners = []
    for assign in itertools.product(range(4), repeat=nt):
        sgn = {t: CHOICES[assign[t]] for t in range(nt)}
        if all(holds(info, sgn) for info in data):
            winners.append(sgn)
    print(f"\nOf all {4**nt} type->sign rules, {len(winners)} give the pentagon for all 120 orderings:")
    for w in winners:
        cflip = [inv[t] for t in range(nt) if w[t][0] == -1]
        sflip = [inv[t] for t in range(nt) if w[t][1] == -1]
        tag = "  <-- paper's B/C rule" if (cflip == [(1, 3, 2, 4)] and sflip == [(1, 2, 4, 3)]) else "  (orientation mirror)"
        print(f"   cos-minus on {cflip},  sin-minus on {sflip}{tag}")

    # Negative tests: the two winners are DISTINCT rules (not one rule relabeled), and the
    # role-swapped version of the paper's rule fails -- so the mirror is not silently identified with
    # the paper's convention.
    def count(rule):
        # rule: dict necklace->(sc,ss); build per-type assignment and count holds/120
        def holds_rule(info):
            def M(name, vec):
                i, j, k, l, c0, s0, ty, rem = info[name]
                sc, ss = rule.get(inv[ty], (1, 1))
                return flip_on_vector(vec, i, j, k, l, sc * c0, ss * s0, rem)
            for t0 in T0_BASIS:
                top = M("ii", M("i", {t0: 1 + 0j}))
                bot = M("v", M("iv", M("iii", {t0: 1 + 0j})))
                if any(abs(top.get(x, 0) - bot.get(x, 0)) > 1e-9 for x in set(top) | set(bot)):
                    return False
            return True
        return sum(holds_rule(info) for info in data)
    paper = {(1, 3, 2, 4): (-1, 1), (1, 2, 4, 3): (1, -1)}
    mirror = {(1, 4, 2, 3): (-1, 1), (1, 3, 4, 2): (1, -1)}
    swapped = {(1, 3, 2, 4): (1, -1), (1, 2, 4, 3): (-1, 1)}   # cos/sin roles swapped on B,C
    assert count(paper) == 120 and count(mirror) == 120 and paper != mirror, "winner identity test"
    assert count(swapped) < 120, "role-swap negative test"
    print(f"\n   negative tests: paper rule {count(paper)}/120, mirror rule {count(mirror)}/120 "
          f"(distinct rules); role-swapped paper rule {count(swapped)}/120 < 120. PASS")


if __name__ == "__main__":
    main()

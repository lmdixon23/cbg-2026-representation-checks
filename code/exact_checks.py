"""
exact_checks.py -- exact-arithmetic re-verification (float/tolerance checks are diagnostics,
not certificates).

Everything here is computed in exact arithmetic with sympy (rationals + symbolic square roots);
no floating point, no tolerance. Establishes C1 (Ptolemy), C3 (sign uniqueness), and C4 (pentagon,
far-commute) in exact finite arithmetic.

Run:  python exact_checks.py    (a few seconds)
"""
import sys, os, itertools
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sympy as sp
from operators import classify, paper_signs, flip_on_vector
from relations import PENT_COORDS, PENT_FLIPS, T0_BASIS


def is_zero(expr):
    """Exact zero test for an algebraic-number expression."""
    e = sp.simplify(sp.expand(expr))
    return e == 0


# ---- C1: Ptolemy identity, symbolic (exact for ALL values, not sampled) ----
def ptolemy_symbolic():
    zi, zj, zk, zl = sp.symbols('zeta_i zeta_j zeta_k zeta_l')
    expr = (zi - zl) * (zj - zk) + (zi - zj) * (zk - zl) - (zi - zk) * (zj - zl)
    ok = sp.expand(expr) == 0
    print(f"C1  Ptolemy z_il z_jk + z_ij z_kl - z_ik z_jl == 0 symbolically: {ok}  (exact, all values)")
    return ok


# ---- exact cos/sin and exact pentagon ----
def cos_sin_exact(vals, i, j, k, l, sc, ss):
    z = lambda a, b: vals[a] - vals[b]
    c = sc * sp.sqrt(z(i, l) * z(j, k) / (z(i, k) * z(j, l)))
    s = ss * sp.sqrt(z(i, j) * z(k, l) / (z(i, k) * z(j, l)))
    return c, s


def pentagon_exact(vals, sign_overrides=None):
    """Exact 5-cycle check on the T0 basis for one zeta-assignment. sign_overrides maps a necklace
    to (sc,ss) to override paper_signs (used by the sign-uniqueness confirmation)."""
    def M(name, vec):
        quad, rem, _ = PENT_FLIPS[name]
        i, j, k, l, neck = classify(quad, vals, PENT_COORDS)
        sc, ss = (sign_overrides.get(neck, (1, 1)) if sign_overrides is not None else paper_signs(neck))
        c, s = cos_sin_exact(vals, i, j, k, l, sc, ss)
        return flip_on_vector(vec, i, j, k, l, c, s, rem)
    for t0 in T0_BASIS:
        top = M("ii", M("i", {t0: sp.Integer(1)}))
        bot = M("v", M("iv", M("iii", {t0: sp.Integer(1)})))
        for x in set(top) | set(bot):
            if not is_zero(top.get(x, 0) - bot.get(x, 0)):
                return False
    return True


def pentagon_all_orderings_exact():
    vals_pool = [sp.Integer(v) for v in (5, 4, 3, 2, 1)]
    ok = 0
    for perm in itertools.permutations(vals_pool):
        vals = {i + 1: perm[i] for i in range(5)}
        ok += pentagon_exact(vals)
    print(f"C4  pentagon (5-cycle) with paper sign rule: {ok}/120 orderings hold EXACTLY")
    return ok == 120


# ---- C3: exact confirmation of the two winning sign rules + refutation of neighbors ----
def sign_uniqueness_exact():
    B, C = (1, 3, 2, 4), (1, 2, 4, 3)          # paper's cos-minus / sin-minus types
    Bm, Cm = (1, 4, 2, 3), (1, 3, 4, 2)        # the orientation-mirror pair
    winners = [{B: (-1, 1), C: (1, -1)}, {Bm: (-1, 1), Cm: (1, -1)}]
    losers = [{B: (-1, 1)},               # only one flip
              {B: (1, -1), C: (-1, 1)},   # swapped roles
              {B: (-1, 1), C: (-1, 1)}]   # both cos-minus
    vals_pool = [sp.Integer(v) for v in (5, 4, 3, 2, 1)]
    def count(ov):
        return sum(pentagon_exact({i + 1: perm[i] for i in range(5)}, ov)
                   for perm in itertools.permutations(vals_pool))
    wins = [count(w) for w in winners]
    losses = [count(l) for l in losers]
    print(f"C3  winning sign rules hold exactly 120/120: {wins}  (both == 120: {all(w==120 for w in wins)})")
    print(f"C3  representative non-winning rules fall short: {losses}  (all < 120: {all(l<120 for l in losses)})")
    return all(w == 120 for w in wins) and all(l < 120 for l in losses)


# ---- C3 exhaustive: EXACT certification that exactly two of 4^6 sign rules win ----
def _float_pent_fails(geo, rule):
    """Float pentagon check used ONLY to locate a candidate failing ordering (cheap)."""
    import cmath
    def M(name, vec):
        i, j, k, l, neck, rem, fvals = geo[name]
        sc, ss = rule.get(neck, (1, 1))
        z = lambda a, b: fvals[a] - fvals[b]
        c = sc * cmath.sqrt((z(i, l) * z(j, k)) / (z(i, k) * z(j, l)))
        s = ss * cmath.sqrt((z(i, j) * z(k, l)) / (z(i, k) * z(j, l)))
        return flip_on_vector(vec, i, j, k, l, c, s, rem)
    for t0 in T0_BASIS:
        top = M("ii", M("i", {t0: 1 + 0j}))
        bot = M("v", M("iv", M("iii", {t0: 1 + 0j})))
        if any(abs(top.get(x, 0) - bot.get(x, 0)) > 1e-9 for x in set(top) | set(bot)):
            return True
    return False


def _exact_pent_fails_on(info, rule):
    """EXACT: does this ordering's pentagon fail under `rule`? (certifying check)."""
    def M(name, vec):
        i, j, k, l, neck, rem, vals = info[name]
        sc, ss = rule.get(neck, (1, 1))
        c, s = cos_sin_exact(vals, i, j, k, l, sc, ss)
        return flip_on_vector(vec, i, j, k, l, c, s, rem)
    for t0 in T0_BASIS:
        top = M("ii", M("i", {t0: sp.Integer(1)}))
        bot = M("v", M("iv", M("iii", {t0: sp.Integer(1)})))
        for x in set(top) | set(bot):
            if not is_zero(top.get(x, 0) - bot.get(x, 0)):
                return True
    return False


def sign_uniqueness_exhaustive_exact():
    """Certify EXACTLY TWO of all 4^6 type->sign rules satisfy the pentagon for all 120 orderings,
    with an EXACT failure witness for every non-winner (no floating tolerance in the certificate).
    Float is used only to *locate* a candidate witness ordering; the verdict is exact."""
    # geometry per ordering, in exact (vals) and float (fvals) form
    exact_geo, float_geo = [], []
    for perm in itertools.permutations((5, 4, 3, 2, 1)):
        ev = {i + 1: sp.Integer(perm[i]) for i in range(5)}
        fv = {i + 1: float(perm[i]) for i in range(5)}
        ge, gf = {}, {}
        for name, (quad, rem, _) in PENT_FLIPS.items():
            i, j, k, l, neck = classify(quad, ev, PENT_COORDS)
            ge[name] = (i, j, k, l, neck, rem, ev)
            gf[name] = (i, j, k, l, neck, rem, fv)
        exact_geo.append(ge); float_geo.append(gf)

    types = [(1, 2, 3, 4), (1, 3, 2, 4), (1, 2, 4, 3), (1, 4, 3, 2), (1, 4, 2, 3), (1, 3, 4, 2)]
    choices = list(itertools.product([1, -1], repeat=2))
    winners, witnesses = [], {}     # witnesses[rule_idx] = ordering index where it fails (float-located)
    for ridx, assign in enumerate(itertools.product(range(4), repeat=6)):
        rule = {types[t]: choices[assign[t]] for t in range(6)}
        wit = None
        for o in range(120):
            if _float_pent_fails(float_geo[o], rule):
                wit = o; break
        if wit is None:
            winners.append(rule)
        else:
            witnesses[ridx] = (rule, wit)
    # EXACT certification
    winners_ok = all(not any(_exact_pent_fails_on(exact_geo[o], w) for o in range(120)) for w in winners)
    # every non-winner: exactly fails its located ordering
    nonwin_ok = all(_exact_pent_fails_on(exact_geo[o], rule) for (rule, o) in witnesses.values())
    n_win, n_non = len(winners), len(witnesses)
    print(f"C3* exhaustive EXACT: {n_win} winners + {n_non} non-winners = {n_win + n_non} (= 4^6 = 4096)")
    print(f"     both winners verified exactly over all 120 orderings: {winners_ok}")
    print(f"     every one of the {n_non} non-winners has an EXACT failure witness: {nonwin_ok}")
    ok = (n_win == 2 and winners_ok and nonwin_ok)
    print(f"     => EXACTLY TWO winning sign rules, certified exactly: {ok}")
    return ok


# ---- C4: far-commutativity, exact on a representative ordering + structural proof note ----
def far_commute_exact():
    vals = {1: sp.Integer(8), 2: sp.Integer(7), 3: sp.Integer(6), 4: sp.Integer(5),
            5: sp.Integer(4), 6: sp.Integer(3), 7: sp.Integer(2), 8: sp.Integer(1), 9: sp.Rational(1, 2)}
    c1, s1 = cos_sin_exact(vals, 1, 2, 3, 4, 1, 1)
    c2, s2 = cos_sin_exact(vals, 5, 6, 7, 8, 1, 1)
    basis = [frozenset({1, 2, 3}), frozenset({1, 3, 4}), frozenset({5, 6, 7}),
             frozenset({5, 7, 8}), frozenset({2, 3, 9})]
    worst_ok = True
    for t in basis:
        a = flip_on_vector(flip_on_vector({t: sp.Integer(1)}, 1, 2, 3, 4, c1, s1, frozenset({1, 3})),
                           5, 6, 7, 8, c2, s2, frozenset({5, 7}))
        b = flip_on_vector(flip_on_vector({t: sp.Integer(1)}, 5, 6, 7, 8, c2, s2, frozenset({5, 7})),
                           1, 2, 3, 4, c1, s1, frozenset({1, 3}))
        for x in set(a) | set(b):
            if not is_zero(a.get(x, 0) - b.get(x, 0)):
                worst_ok = False
    print(f"C4  far-commutativity exact on representative ordering: {worst_ok}")
    print("    structural proof: flips in vertex-disjoint quadrilaterals act on disjoint triangle-")
    print("    basis vectors (supp {ijk,ikl,ijl,jkl} disjoint), so as operators they live on")
    print("    complementary coordinate blocks and commute; the exact check confirms one instance.")
    return worst_ok


def main():
    print("== Exact-arithmetic re-verification (sympy; no floating point, no tolerance) ==\n")
    r = [ptolemy_symbolic(), pentagon_all_orderings_exact(), sign_uniqueness_exact(),
         sign_uniqueness_exhaustive_exact(), far_commute_exact()]
    print(f"\nAll exact checks passed: {all(r)}")


if __name__ == "__main__":
    main()

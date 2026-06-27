"""
check_flip_sequence.py -- exact/combinatorial certificate for the Delaunay flip sequences that
exact_kernel.py applies (Checker B of verification/CHECKER-DESIGN-C8prime-C9-v1.md; partial discharge
of the R2 obligation: the C9 kernel result rests on a float-located flip sequence).

This script does NOT trust the floating geometry's verdict. For each inner pure-braid generator it:
  (B1) takes the flip sequence the pipeline detects, then INDEPENDENTLY certifies, by exact
       combinatorial predicates (frozensets only; no float, no tolerance), that the sequence is a legal
       2->2 flip path -- every step removes a diagonal whose two triangles are BOTH present, replaces it
       with the other diagonal -- and that it returns to the initial triangulation (a closed loop);
  (B2-lite) re-detects the sequence at several sampling resolutions and asserts the combinatorial
       sequence is identical -- robustness to the float sampling, ruling out a near-degenerate
       in-circle mis-detection.

What remains OPEN (CHECKER-DESIGN-..., Checker B2-full): a purely combinatorial derivation of each
generator's flip word from the braid word with NO geometry at all. Until that exists C9 stays
CONDITIONALLY CERTIFIED; this script upgrades the condition from "unverified float sequence" to
"verified legal, closed, and resolution-robust float sequence".

Run:  python check_flip_sequence.py
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from exact_kernel import flip_sequence
from float_kernel import gen_words


def certify_legal_loop(init, seq):
    """EXACT combinatorial: is `seq` a legal 2->2 flip path from `init` back to `init`?
    Each flip names quad (i,j,k,l) and the removed diagonal; the two triangles sharing that diagonal
    must both be present, and are replaced by the two triangles of the other diagonal. Hardened
    (2026-06-27, hardening pass) toward incapable-of-falsely-reporting form: the two inserted
    triangles must be genuinely new (not already in the triangulation), and the triangle count must be
    invariant at every step (a legal 2->2 flip removes exactly two and adds exactly two)."""
    tris = set(init)
    n0 = len(tris)
    for step, (i, j, k, l, removed) in enumerate(seq):
        ijk, ikl = frozenset({i, j, k}), frozenset({i, k, l})
        ijl, jkl = frozenset({i, j, l}), frozenset({j, k, l})
        if set(removed) == {i, k}:                 # forward flip ik -> jl
            old, new = (ijk, ikl), (ijl, jkl)
        elif set(removed) == {j, l}:               # inverse flip jl -> ik
            old, new = (ijl, jkl), (ijk, ikl)
        else:
            return False, f"step {step}: removed {set(removed)} is not a diagonal of quad {{{i},{j},{k},{l}}}"
        if old[0] not in tris or old[1] not in tris:
            return False, f"step {step}: illegal flip -- removing {set(removed)} but its two triangles are not both present"
        tris.discard(old[0]); tris.discard(old[1])
        if new[0] in tris or new[1] in tris:       # a 2->2 flip must introduce two GENUINELY NEW triangles
            return False, f"step {step}: illegal flip -- an inserted triangle already exists (not a clean 2->2 flip)"
        tris.add(new[0]); tris.add(new[1])
        if len(tris) != n0:                        # triangle count is invariant under a legal 2->2 flip
            return False, f"step {step}: triangle-count changed to {len(tris)} (expected {n0})"
    return (tris == set(init)), ("closed loop" if tris == set(init) else "does NOT return to the initial triangulation")


def seq_key(seq):
    return tuple((i, j, k, l, tuple(sorted(removed))) for (i, j, k, l, removed) in seq)


def check(n):
    print(f"[n={n}] certifying {len(gen_words(n))} inner-generator flip sequences:")
    ok_all = True
    for ij, word in gen_words(n).items():
        init, final, seq = flip_sequence(n, word, steps=600)
        ok_loop, why = certify_legal_loop(init, seq)
        # B2-lite: same combinatorial sequence across resolutions (rules out near-degenerate mis-detection)
        keys = {seq_key(flip_sequence(n, word, steps=s)[2]) for s in (300, 600, 1200)}
        robust = (len(keys) == 1)
        ok = ok_loop and robust
        ok_all = ok_all and ok
        print(f"   A_{ij}: {len(seq):2d} flips -> legal closed loop: {ok_loop} ({why}); "
              f"resolution-robust 300/600/1200: {robust}  {'OK' if ok else 'X'}")
        assert ok_loop, f"A_{ij}: flip sequence is not a legal closed loop ({why})"
        assert robust, f"A_{ij}: flip sequence is NOT resolution-robust (float sampling artifact)"
    return ok_all


def main():
    print("== Exact/combinatorial flip-sequence certificate (B1 legal closed loop + B2-lite robustness) ==\n")
    check(3)
    print()
    check(4)
    print("\n=> every inner-generator flip sequence (n=3,4) is an exact-combinatorially LEGAL closed loop and is")
    print("   robust to sampling resolution. C9's float-sequence condition is thereby VERIFIED (legal + robust);")
    print("   eliminating the geometry entirely (Checker B2-full: combinatorial flip word from the braid) is OPEN.")


if __name__ == "__main__":
    main()

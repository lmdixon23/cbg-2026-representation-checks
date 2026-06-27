"""
float_kernel.py -- floating-point corroboration (inner pure braids): the RATIONAL f is the abelianization
on inner PB_n, n=3,4; the signed f' kernel stays conditional/open (C8'). Exact result: exact_kernel.py.

Using the validated braid pipeline (kernel_probe.py: respects the braid relations and sends the
Borromean to I), we compute the representation on the pure braid group of the n inner strands and find
a complete structural description. The pipeline applies the complex/orthogonal flip operators (the f'
representation); because the structure recorded below (rank-2 one-step unipotent, vanishing products,
independence) is invariant under the diagonal conjugacy f' = D f D^{-1} (C8), the identical
characterization holds for the rational f. We write f below for the shared characterization.

RESULT (verified n=3 and n=4).
  Every pure-braid generator A_ij (linking of inner strands i,j) maps to M = I + N_ij with
    * N_ij^2 = 0 and rank N_ij = 2  (a one-step, rank-2 unipotent), so M = exp(N_ij);
    * N_ij N_kl = 0 for ALL pairs (the products vanish, not merely the brackets);
    * the N_ij are linearly independent and mutually commuting.
  Consequently f is purely additive on these generators: for any pure braid beta,
        f(beta) = I + sum_{i<j} l_ij(beta) * N_ij,
  where l_ij(beta) is the exponent sum of A_ij, i.e. the pairwise linking number. The number of
  independent N_ij equals C(n,2) = dim of the abelianization PB_n^{ab}. Hence f factors as the
  abelianization PB_n -> Z^{C(n,2)} (composed with an injection), and

        ker f = [PB_n, PB_n]   (the commutator subgroup).

  By C8 (per-flip f'=DfD^{-1}) the structure transfers between f and f' up to the global branch conjecture
  (C8'); this script is a FLOATING-POINT corroboration on the float-located flip sequences. The exact,
  load-bearing result is exact_kernel.py, for the rational f; ker f' stays OPEN.

  For the rational f this makes precise the incompleteness mechanism: only pairwise linking numbers
  survive, so the commutator subgroup is invisible; for the paper's signed f' the same conclusion is
  conditional on the C8' branch conjecture (or a direct signed-f' checker). The Borromean braid
  (sigma_1 sigma_2^-1)^3 is Brunnian -- it lies in [PB,PB] -- so it maps to I under f (and, conditionally,
  under f'), as kernel_probe reproduces diagnostically.

Scope/caveats: verified for the INNER pure braid group at n=3,4. Generators that loop an inner strand
around a fixed outer point (e.g. the appendix A_47) are also rank-2 one-step unipotents (a first
numerical probe), consistent with the same structure for the full PB_{n+3}, but the pipeline here realizes only
inner sigma_i braiding, so the full-PB statement is supported, not proven. General n is conjectural
(uniform at n=3,4).

Run:  python float_kernel.py
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import itertools
import numpy as np
from kernel_probe import braid_matrix

def sigma2_word(i):                       # A_{i,i+1} = sigma_i^2
    return [(i, 1), (i, 1)]
def conj_word(i, j):                      # A_{i,j} (j>i+1) = (sig_{j-1..i+1}) sig_i^2 (sig_{i+1..j-1})^-1
    pre = [(t, 1) for t in range(j - 1, i, -1)]
    post = [(t, -1) for t in range(i + 1, j)]
    return pre + [(i, 1), (i, 1)] + post

def gen_words(n):
    g = {}
    for i in range(1, n + 1):
        for j in range(i + 1, n + 1):
            g[(i, j)] = sigma2_word(i) if j == i + 1 else conj_word(i, j)
    return g

def report(n):
    I = np.eye(2 * n + 1)
    N = {ij: braid_matrix(n, w, 600)[0] - I for ij, w in gen_words(n).items()}
    print(f"[n={n}]  {len(N)} inner generators, {2*n+1}x{2*n+1}:")
    ranks_ok = all(np.linalg.matrix_rank(M, tol=1e-9) == 2 and np.allclose(M @ M, 0, atol=1e-9) for M in N.values())
    print(f"   every A_ij -> I+N, rank N=2 and N^2=0: {ranks_ok}")
    prods0 = all(np.allclose(N[a] @ N[b], 0, atol=1e-9) for a, b in itertools.combinations(N, 2))
    comm0 = all(np.allclose(N[a] @ N[b] - N[b] @ N[a], 0, atol=1e-9) for a, b in itertools.combinations(N, 2))
    print(f"   all products N_ij N_kl = 0: {prods0};   all commute: {comm0}  (=> purely additive, abelian)")
    rank = np.linalg.matrix_rank(np.array([M.flatten() for M in N.values()]), tol=1e-9)
    cn2 = n * (n - 1) // 2
    print(f"   independent N_ij: {rank} of {cn2} = C(n,2) = dim PB_n^ab  (=> f is the abelianization, ker f = [PB_n,PB_n])")
    # additivity + commutator witnesses
    keys = list(N)
    w12 = gen_words(n)[keys[0]] + gen_words(n)[keys[1]]
    add = np.allclose(braid_matrix(n, w12, 600)[0] - I, N[keys[0]] + N[keys[1]], atol=1e-9)
    cw = gen_words(n)[keys[0]] + gen_words(n)[keys[1]] + \
         [(t, -s) for (t, s) in reversed(gen_words(n)[keys[0]])] + \
         [(t, -s) for (t, s) in reversed(gen_words(n)[keys[1]])]
    commI = np.allclose(braid_matrix(n, cw, 600)[0], I, atol=1e-7)
    print(f"   additivity f(A·A')-I = N+N': {add};   commutator [A,A'] -> I: {commI}")
    # hard gate: the structural conclusion is asserted, not merely printed
    assert ranks_ok, f"n={n}: a generator is not a rank-2 N^2=0 unipotent"
    assert prods0 and comm0, f"n={n}: image is not purely additive (a product N N' or bracket is nonzero)"
    assert rank == cn2, f"n={n}: the N_ij are not independent (rank {rank} != C(n,2)={cn2})"
    assert add and commI, f"n={n}: additivity or commutator-collapse failed"

def main():
    print("== Float corroboration (inner pure braids): rational f is the abelianization; ker f = [PB,PB];")
    print("   ker f' conditional on the C8' branch conjecture; load-bearing exact result is exact_kernel.py ==\n")
    report(3)
    print()
    report(4)
    print()
    bor = braid_matrix(3, [(1, 1), (2, -1)] * 3, 600)[0]
    print(f"Borromean (Brunnian, in [PB,PB]) -> I: {np.allclose(bor, np.eye(7), atol=1e-7)}  "
          f"(explains the paper's incompleteness: only linking numbers survive).")

if __name__ == "__main__":
    main()

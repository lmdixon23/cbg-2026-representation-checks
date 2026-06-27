# A verification note on the complex-orthogonal representation of the colored braid groupoid

Status: draft
Code and verification ledgers for the paper (arXiv link upon posting).
Prompted by and commenting on arXiv:2606.20473 (Rohozhkin, *Invariants of the Colored Braid Groupoid*).

This is a verification-and-corrections study of the rational representation `f` and the
complex-orthogonal representation `f'` of the colored braid groupoid, with three original results:
the per-flip diagonal similarity `f' = D f D^{-1}` (C8; the global transfer to the paper's B/C-signed
`f'` is conjectural — C8′), the identification of the **rational** invariant with the abelianization /
pairwise linking numbers so that `ker f = [PB,PB]` on the inner pure braid group for `n = 3,4` (C9, in
exact arithmetic given the Delaunay-located flip sequence; `ker f'` open), and an erratum in the
paper's two displayed worked-example matrices (C10).

## Reproduce everything

    pip install -r requirements.txt
    bash run_all.sh          # exact certificate first; ~2 min total (exact_kernel.py dominates)

On Windows without WSL, run the scripts with the same interpreter your packages are installed in:

    pip install -r requirements.txt
    @("code\exact_checks.py","code\diagonal_form.py","code\pentagon_typo.py",
      "code\appendix_reconstruction.py","code\appendix_full_matrix.py","code\relations.py",
      "code\sign_uniqueness.py","code\kernel_probe.py","code\float_kernel.py",
      "code\exact_kernel.py") | % { Write-Host "`n### $_"; python $_ }

Run scripts from the repository root. The finite claims are verified in exact arithmetic
(`sympy`); the floating-point scripts are fast smoke tests / reconstructions that the exact
checks supersede.

## Claim-to-script map

| Claim | Script | What to look for |
|---|---|---|
| C1: `cos²φ+sin²φ=1 ⇔` Ptolemy identity `z_il z_jk + z_ij z_kl = z_ik z_jl`; image in `SO(2n+1,ℂ)` | `code/exact_checks.py`; `code/operators.py` | "C1 Ptolemy ... == 0 symbolically: True"; single-block complex-orthogonal True |
| C2: main-text pentagon typo (3rd RHS factor `φ₂₃₄₅`→`φ₁₂₄₅`) | `code/pentagon_typo.py` | printed form ‖·‖≈0.26 FAILS; appendix form ≈1e-16 HOLDS |
| C3: sign rule essentially unique (exactly 2 of `4⁶`) | `code/exact_checks.py`; `code/sign_uniqueness.py` | "EXACTLY TWO winning sign rules, certified exactly: True"; mirror/role-swap negatives |
| C4: four groupoid relations (pentagon 120/120, far-commute) | `code/exact_checks.py`; `code/relations.py` | pentagon 120/120 exact; far-commute 0 over all 8! orderings |
| C5: two printed flip matrices reproduce (convention lock) | `code/appendix_reconstruction.py` | `A_4615` cos 0.79056942, sin 0.61237244; `A_1524` cos 1.06066017, sin 0.35355339·i |
| C6: rational `f` unit-diagonal/det-1 on `b₄,₇` | `code/appendix_reconstruction.py` | reversed order: max\|diag−1\|≈4e-16, det 1 |
| C8: **per-flip** diagonal similarity `f'=DfD⁻¹` via `fᵀQf=Q` (global signed-`f'` transfer conjectural — C8′) | `code/diagonal_form.py` | "f^T Q_p f − Q_t = 0 -> True"; `b₄,₇` Q-orthonormalized f' orthogonal/unit-diag/det 1 |
| C9: **rational** invariant = abelianization, `ker f=[PB,PB]` (inner `n=3,4`, given the float-located flip sequence; `ker f'` open) | `code/exact_kernel.py`; `code/float_kernel.py` | n=3,4: rank-2 N²=0, all products 0, independent `C(n,2)`; Borromean → I |
| C10: appendix displayed result matrices erroneous | `code/appendix_full_matrix.py` | displayed det 617.52 / 49 (I+skew); factor products orthogonal, det 1, rank-2 unipotent |
| Braid → Delaunay-flip pipeline (diagnostic tool) | `code/kernel_probe.py` | checkpoints (braid relation, far-commutation, Borromean) all → I |

## Verification methodology

The finite claims (C1, C3, C4) and the structural results (C8, C9 at inner `n=3,4`) are checked in
exact arithmetic with `sympy` — no floating point, no tolerance. `verification/STATUS-LEDGER.md`
records the status and evidence type of every claim. The results were cross-checked by an independent
re-derivation and by a second, decorrelated verification pass.

## Independent author check

`BLIND_CHECK.py` recomputes the headline quantities — the two printed flip matrices, the
Ptolemy/orthogonality identity, the single-flip `fᵀQf=Q` form-invariance, and the `b₄,₇`
determinant/unit-diagonal/rank-2 properties — from the formulas alone, importing nothing from `code/`.
It is a second, independent implementation, so a transcription or convention error in `code/` would
surface as a disagreement. Spec and completion record: `verification/BLIND_CHECK.md`.

    python BLIND_CHECK.py     # exit 0 and "BLIND CHECK PASSED"

## Scope and caveats

C9 is established exactly for the **inner** pure braid group at `n = 3, 4`; the general-`n` and full
`PB_{n+3}` statements are conjectural. C2 and C10 are claims about the paper's printed text; the source
(arXiv:2606.20473v1) has been read directly and independently re-read in a fresh context, confirming both,
so they are CERTIFIED (C5's flip-matrix values are independently re-derived from the operator formulas and
likewise do not depend on transcription).

## Cite

See `CITATION.cff`.

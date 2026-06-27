# STATUS-LEDGER — cbg-2026-representation-checks verification

Single source of truth for the status of every load-bearing claim in the accompanying note.
Finite claims are verified in exact arithmetic (`sympy`); floating-point scripts are fast
reconstructions the exact checks supersede. Status vocabulary: CERTIFIED · CONDITIONALLY CERTIFIED ·
DIAGNOSTIC ONLY · OPEN. Evidence type: algebraic identity · finite exhaustive (exact) ·
checker-verified · structural proof · strategic inference.

A claim is CERTIFIED when an exact or executed check establishes it, the check has been reproduced by
an independent, decorrelated review, and it carries no unmet caveat. CONDITIONALLY CERTIFIED marks a
claim whose computation is verified but which rests on an external dependency (transcription of the
source paper) or holds only over a declared sub-domain.

**Status (2026-06-27).** An independent, decorrelated review of the artifact adopted two Serious
soundness findings — the global `f'=DfD⁻¹` transfer to the paper's B/C-signed operator was overstated
(R1), and the kernel exactness rests on a float-located flip sequence (R2) — plus two Moderate (R3
pipeline wording, R4 the constant `≈617.52`); C8 and C9 were scoped accordingly (below), and a re-review
of the corrected source then closed with no remaining soundness finding. Separately, the C2/C10
transcription caveat is discharged: arXiv:2606.20473v1 was read directly and independently re-read, both
confirming the main-text pentagon typo (eq 2.8 third RHS factor `φ₂₃₄₅`; appendix `φ₁₂₄₅`) and that the
two displayed result matrices are unit-diagonal `I+skew` (not orthogonal). C2 and C10 are CERTIFIED.
Full adjudication history is kept in `_local/STATUS-LEDGER-history.md`.

| # | Claim | Evidence type | Status | Script |
|---|---|---|---|---|
| C1 | `cos²φ+sin²φ = 1 ⇔` Ptolemy/Plücker identity `z_il z_jk + z_ij z_kl = z_ik z_jl`; each flip operator is complex-orthogonal with determinant `+1`, so the image lies in `SO(2n+1,ℂ)` | algebraic identity | CERTIFIED | `code/exact_checks.py`, `code/operators.py` |
| C2 | the main-text pentagon equation has a typo: its 3rd right-hand factor `φ₂₃₄₅` should be `φ₁₂₄₅` (the appendix is correct) | checker-verified + structural + source-verified | CERTIFIED | `code/pentagon_typo.py` |
| C3 | the quadrilateral-type sign rule is essentially unique: of `4⁶` rules exactly two satisfy the pentagon for all 120 orderings — the paper's B/C rule and its orientation mirror | finite exhaustive (exact) | CERTIFIED | `code/exact_checks.py`, `code/sign_uniqueness.py` |
| C4 | all four groupoid relations hold (identity, 4-cycle, pentagon 120/120, far-commutativity) | finite exhaustive (exact) + structural proof | CERTIFIED | `code/exact_checks.py`, `code/relations.py` |
| C5 | the two printed appendix flip matrices reproduce exactly under the positional-subscript / principal-sqrt convention | checker vs source + formula re-derivation | CERTIFIED | `code/appendix_reconstruction.py` |
| C6 | the rational representation `f` has unit diagonal and determinant 1 on the worked pure braid `b₄,₇` | checker-verified | CERTIFIED | `code/appendix_reconstruction.py` |
| C7 | unit diagonal ⇒ constant trace `2n+1` ⇒ the naïve trace is not a knot/link invariant (the paper's stated obstruction) | strategic inference | DIAGNOSTIC ONLY | (note remark) |
| C8 | `f` preserves the diagonal form `Q(Δ_pqr)=(ζ_pq ζ_pr ζ_qr)^{-1}` (`fᵀ Q f = Q`, exact), and dividing each rational flip by `D=diag(Q^{1/2})` yields the complex-orthogonal flip, so `f'=DfD⁻¹` **per flip** | algebraic identity (exact) | **CERTIFIED (per-flip)** | `code/diagonal_form.py` |
| C8′ | the paper's **global B/C sign rule** realizes this `Q^{1/2}` branch around every closed loop, so unit-diagonal / det-1 / kernel transfer from `f` to the paper's **own** signed `f'` | strategic inference (branch bookkeeping) | **OPEN DESIGN** (explicit conjecture; example-confirmed on `b₄,₇` only) | `code/diagonal_form.py` |
| C9 | the **rational** invariant `f` is the abelianization: each `A_ij ↦ I+N_ij` with `N²=0`, rank 2, all products `N_ij N_kl=0`, the `N_ij` independent ⇒ `f(β)=I+Σ ℓ_ij(β) N_ij`, so `ker f=[PB,PB]` | additive structure, exact (sympy) **given the float-located flip sequence** | **CONDITIONALLY CERTIFIED** (inner `PB_n`, `n=3,4`; exact arithmetic conditional on the Delaunay-located flip sequence — now verified legal + resolution-robust by Checker B (`code/check_flip_sequence.py`, ran clean 2026-06-27; all 9 inner generators, n=3,4), still geometry-located pending B2-full). `ker f'=[PB,PB]` is **OPEN IMPLEMENTATION** (needs the C8′ bridge or a direct exact `f'` check); general `n` / full `PB_{n+3}` OPEN | `code/exact_kernel.py`, `code/float_kernel.py` |
| C10 | the paper's two displayed appendix result matrices are erroneous (unit-diagonal but `I + skew`, determinants ≈617.52 and 49, not orthogonal); the products of its printed factors are correct (orthogonal, det 1, rank-2 unipotent) | full-matrix audit + source-verified | CERTIFIED | `code/appendix_full_matrix.py` |
| O1 | complex `f'` product is unit-diagonal on pure braids | algebraic, via C8 | RESOLVED for the `Q`-orthonormalized `f'`; for the paper's signed `f'` conditional on C8′ | `code/diagonal_form.py` |
| O2 | the Borromean braid `(σ₁σ₂⁻¹)³` maps to `I` | pipeline + structural (Brunnian ∈ `[PB,PB]`) | RESOLVED | `code/kernel_probe.py` |
| O3 | the kernel of `f'` on the inner pure braid group | C9 | `ker f` characterized (inner `n=3,4`, conditional-exact per C9); `ker f'` OPEN (C8′/C9) | `code/exact_kernel.py` |

## Open obligations
- **C8′ (R1):** an exact branch checker for the B/C sign rule — or a direct global proof that the
  paper's signed `f'` equals `DfD⁻¹` around every closed loop — to lift the global `f'` transfer from
  conjecture to certified.
- **C9 (R2):** `code/check_flip_sequence.py` (Checker B) provides **B1** — an exact combinatorial certificate
  that each inner generator's detected flip sequence is a legal closed flip loop — plus **B2-lite**
  (resolution-robustness, ruling out near-degenerate mis-detection). **Ran clean (2026-06-27):** all 9
  inner-generator sequences (n=3,4) are exact-combinatorially legal closed loops and resolution-robust at
  300/600/1200 samples, so C9's condition is now *verified legal + robust* (no longer unverified). **Still open:** B2-full — a purely
  combinatorial flip word derived from the braid with no geometry — to lift C9 to fully exact
  (`verification/CHECKER-DESIGN-C8prime-C9-v1.md`).
- **C9 / signed `f'` (R3):** a direct exact computation of the **signed** `f'` kernel on inner
  `PB_{3,4}`, rather than transfer through the C8′ conjecture.
- C9 for general `n` / the full `PB_{n+3}`: requires a uniform proof that `N_ij N_kl = 0` for all
  pairs (the products, not merely the brackets, vanish).
- C2 / C10 — **DISCHARGED (2026-06-27).** The source (arXiv:2606.20473v1) was read directly and independently
  re-read; both confirm the pentagon typo (eq 2.8 third RHS factor `φ₂₃₄₅`, appendix `φ₁₂₄₅`) and the
  `I+skew` result matrices, so C2/C10 are now CERTIFIED. (A human glance at eq 2.8 and the two result
  matrices is the trivial belt-and-suspenders.)

> The checkers that would discharge the C8′ (R1) and C9 (R2) obligations — an exact B/C-branch
> certificate and an exact flip-sequence certificate — are scoped at Gate A in
> `verification/CHECKER-DESIGN-C8prime-C9-v1.md` (design only; not yet implemented).

## Reproducibility
Every row is reproduced by the named script; `run_all.sh` runs them in claim order. Environment:
`numpy`, `scipy`, `sympy` (see `requirements.txt`). The exact checks use no floating point and no
tolerance.

## Changelog
- v2.2 (2026-06-27) — an independent, decorrelated review of v0.4 returned CLOSED; it reproduced
  `appendix_full_matrix.py` + `check_flip_sequence.py`, confirming C2/C10 and all 9 Checker-B sequences.
  Non-blocking edits applied (a citation split; `certify_legal_loop` hardened with new-triangle and
  cardinality-invariant guards, re-ran 9/9 clean). No status row changes.
- v2.1 (2026-06-27) — C2/C10 transcription caveat discharged → both promoted to CERTIFIED. The source was
  read directly and independently re-read; C10's verdict is also a model-independent structural fact (a
  product of complex-orthogonal factors is complex-orthogonal, but `I+`nonzero-real-skew is not).
- v2 (2026-06-26) — reopened by a decorrelated review; adopted R1 (→ new row C8′, OPEN DESIGN), R2 (C9
  conditional-exact on the float-located flip sequence; `ker f'` OPEN IMPLEMENTATION), R3, R4. A re-review
  of the corrected source returned CLOSED after propagation fixes.
- v0.3 src-rev c (2026-06-26) — additive: added `code/check_flip_sequence.py` (Checker B) toward the C9/R2
  obligation and wired it into `run_all.sh`. Ran clean 2026-06-27 (all 9 inner generators legal closed loops
  and resolution-robust); C9's condition now verified legal + robust, B2-full still open.
- v1 (2026-06-24) — created.

Full process detail (review rounds, raw output) is kept privately in `_local/STATUS-LEDGER-history.md`.

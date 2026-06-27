# CHECKER DESIGN v1 (2026-06-26) — discharging the two reopened obligations (C8′, C9)

Status: **DESIGN** (Gate-A; not yet implemented). Track: math-proof. Governs the two checkers whose absence
kept the 2026-06-26 review open (R1 → C8′, R2 → C9). Each section gives a contract, the approach, the
validation checkpoints, the failure modes, and the evidence upgrade it would earn. Implementation is a future
Gate-B/C cycle.

The two obligations are independent and can be built in either order. **Checker B (C9) is the higher-leverage
and lower-risk one** — it removes a floating-point dependency from an otherwise-exact result; build it first.

---

## Checker A — exact B/C-branch certificate (lifts C8′ from OPEN DESIGN)

### Contract
- **Claim under test:** the paper's **B/C-signed** operator `f'_BC` equals the `Q`-orthonormalized conjugate
  `D f D^{-1}` (with a single fixed `D=diag(Q_{t_0}^{1/2})`) on the pure braid group — i.e. the B/C sign
  rule is exactly the `Q^{1/2}` branch bookkeeping, so unit-diagonal / det-1 / kernel transfer from `f` to
  the paper's own `f'`.
- **Soundness rule (incapable-of-falsely-reporting form):** the checker reports PASS for a generator `A` only
  if `f'_BC(A) - D f(A) D^{-1}` simplifies to the exact zero matrix in `sympy` over `ℚ(√·)`; any nonzero
  entry → FAIL with that entry as witness.
- **Scope:** the inner pure-braid generators at `n=3,4` (finite, exact). A PASS gives **CONDITIONALLY
  CERTIFIED (inner `PB_{3,4}`)** for C8′; the **general-`n` / full-groupoid** statement needs the cocycle
  argument below and stays OPEN.

### Approach
1. **Per-flip branch extraction (exact).** For each flip, with `ζ_l=l` so the `Q` values and `cos/sin`
   radicands are rationals, compute (i) the rational block `f_flip`, (ii) the B/C-signed block
   `[[c,s],[-s,c]]` with the paper's `(-1)^s,(-1)^r` signs, (iii) `D_p f_flip D_t^{-1}` with
   `D=diag(Q^{1/2})` on the touched triangles. The single-flip identity already holds **up to** a per-triangle
   sign `σ(t)∈{±1}`; solve for the sign relation the B/C rule imposes among the flip's four triangles.
2. **Direct loop check (finite, exact) — the certificate.** For each pure-braid generator `A_ij` (the worked
   `b₄,₇=A₄₇`, `b₅,₆=A₅₆`, and the inner `PB_{3,4}` generators), build `f'_BC(A)` from the exact signed flip
   blocks along the (certified, Checker-B) flip sequence, and `D f(A) D^{-1}` with one fixed `D_{t_0}`, and
   assert `sympy.simplify(f'_BC(A) - D f(A) D^{-1}) == 0`. PASS on all generators ⇒ `f'_BC = D f D^{-1}` on
   the tested `PB`.
3. **Cocycle framing (for the general statement, design-only here).** The per-flip signs define a branch
   1-cochain over the flip groupoid; global consistency = that cochain is a **coboundary**, i.e. its product
   around each defining relation (pentagon 5-cycle, far-commutation) is `+1`. Reducing C8′ to a finite check
   on the relations (exact, like C3/C4) is the route to a general-`n` proof; specify it but defer.

### Validation checkpoints
- On `b₄,₇`: `diagonal_form.py` already recovers `Q` numerically and gets `f'(b₄,₇)` orthogonal/unit-diag/det-1
  (`~1e-16`). Checker A must reproduce this **exactly** (`f'_BC(b₄,₇) = D f(b₄,₇) D^{-1}` symbolically).
- Negative test: perturb one B/C sign (use the role-swapped rule) and confirm the equality **fails** — the
  checker must distinguish the correct branch from a wrong one.

### Failure modes
- **F-A1 (the conjecture is false).** `f'_BC ≠ D f D^{-1}` for some generator → C8′ is refuted, not certified;
  the kernel transfer to the signed `f'` genuinely fails and the paper keeps `ker f'` open. (This is a real
  possible outcome; the checker is built to be able to report it.)
- **F-A2 (radical simplification).** `sympy.simplify` may not reduce nested `√(rational)` combinations to 0
  automatically; mitigate with `nsimplify`/`radsimp`/`sqrtdenest`, or square-and-compare with sign tracking.
- **F-A3 (`ζ_l=l` is non-generic).** A PASS at `ζ_l=l` is weaker than symbolic `ζ`; strengthen by repeating
  with a second integer assignment and, if tractable, symbolic `ζ`.

---

## Checker B — exact flip-sequence certificate (lifts C9 from "conditional exact" to exact)

### Contract
- **Claim under test:** the flip sequence that `exact_kernel.py` applies (currently **located by floating
  Delaunay geometry**) is (i) a legal `2→2` flip path that returns to the initial triangulation, and (ii)
  the path of the **intended braid** `β`. With both, `f(β)` is the exact operator of `β` and the C9
  structure (`A_ij↦I+N_ij`, `N²=0`, products 0, independent) is exact, not conditional.
- **Soundness rule:** report PASS only if every step passes an **exact combinatorial/geometric predicate**
  (no float, no tolerance) and the sequence equals an **independently derived** combinatorial flip word for
  `β`; any mismatch → FAIL with the offending step.

### Approach (two components; B2 is the gold standard)
1. **B1 — path-legality checker (exact combinatorial).** Emit the flip sequence as data (`list of
   (quad, removed_diagonal, added_diagonal)`). Walk it from the initial triangulation, asserting at each step:
   the removed diagonal is an edge of the current triangulation; the four vertices form the quad; the added
   diagonal is the other diagonal; update the triangulation. Assert it returns to the initial triangulation
   (closed loop). This is pure set/combinatorics — no geometry, no float.
2. **B2 — combinatorial braid realization (removes float entirely).** Derive the flip sequence of each braid
   generator `σ_i` (a position swap of inner strands `i,i+1`) **combinatorially** from the braid word and the
   fixed initial layout — the local crossing produces a determined flip pattern — independent of the float
   simulation. Feed **that** sequence to the exact rational operators. The result is fully exact. Cross-check:
   the combinatorial sequence must match `exact_kernel.py`'s float-detected one on every tested case.
3. **B3 — exact-predicate fallback (validation, not the certificate).** If B2's combinatorial derivation is
   deferred, replace the float in-circle/orientation tests with **exact** predicates (determinant signs) at
   rational sample points along the motion, to verify the detected sequence; this narrows but does not fully
   remove the geometric dependency (sampling between flips), so it yields CONDITIONALLY CERTIFIED, while B2
   yields CERTIFIED.

### Validation checkpoints
- The path-legality walker (B1) must reproduce the initial triangulation for every checkpoint braid
  (`σσ⁻¹`, `σ²σ⁻²`, braid relation, far-commute, Borromean, and the `PB_{3,4}` generators).
- B2's combinatorial sequence must equal the float sequence on those cases (exact agreement); a mismatch is a
  finding (the float detection missed or added a flip).
- Re-running `exact_kernel.py` on the B2 sequence must give the **same** `N_ij` (so the C9 conclusions are
  unchanged) — now with no float in the provenance.

### Failure modes
- **F-B1 (float mis-detection).** A near-degenerate in-circle event made the float pipeline miss/insert a
  flip, so the "exact" matrix was of the wrong braid. B2 (or B3) catches it; this is exactly the risk R2
  named.
- **F-B2 (combinatorial derivation is hard).** The `σ_i` flip pattern may be intricate for inner-outer
  generators; scope B2 to the inner `σ_i` first (where C9 is claimed), defer inner-outer.
- **F-B3 (path-dependence worry — resolved).** `f` is a groupoid representation, so any two legal flip paths
  between the same triangulations give the **same** operator; thus B1+B2 (a legal path that is `β`'s path) is
  sufficient — the operator does not depend on which legal realization is used.

---

## Evidence upgrades if both pass
- C8′: OPEN DESIGN → **CONDITIONALLY CERTIFIED (inner `PB_{3,4}`)**; with the cocycle relation-check, a route
  to general `n`.
- C9: CONDITIONALLY CERTIFIED (float-located) → **CERTIFIED (inner `PB_{3,4}`)** for the rational `ker f`; and,
  combined with a passing Checker A, **CERTIFIED** for the signed `ker f'` on the same domain.
- Re-promotion still routes through the workflow: implement (Gate B) → run (Gate C) → decorrelated
  review/closure round (`evidence-standard-v1.md`) → ledger promotion. This design does not itself promote anything.

## Changelog
2026-06-26 v1 — created at Gate A after the decorrelated review adopted R1 (C8′) and R2 (C9). Designs only.
2026-06-26 — **Checker B implemented** as `code/check_flip_sequence.py`: B1 (exact combinatorial legal-closed-
loop certificate) + B2-lite (resolution-robustness). Gate B done; pending its own Gate-C run + a decorrelated
review before it can lift C9's condition in the ledger. **B2-full** (combinatorial flip word from the braid
word, zero geometry) and **Checker A** (C8′ branch certificate) remain design-only.
2026-06-27 — **Checker B Gate-C run: clean.** All 9 inner-generator flip sequences (n=3,4) certified as legal
closed loops and resolution-robust (300/600/1200 samples). C9's float-sequence condition recorded *verified
(legal + robust)* in the ledger; C9 remains CONDITIONALLY CERTIFIED pending B2-full (geometry-free flip word)
and a decorrelated review of Checker B itself. B2-full and Checker A remain design-only.

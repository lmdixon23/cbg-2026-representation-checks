# Evidence standard v1 (2026-06-12)

Repo-binding copy from the lab standard
(00-lab/methodology/evidence-standard-v1.md). On copying the template:
fill the Instantiation section as artifacts appear; this file ships
public at flip.

## Principle

Evidence proportional to load. An artifact's tier is set by what a
defect in it would falsify, never by its size, age, or author.

## Tiers

- T1, certifying: a false output makes a published claim false.
  Checkers, certifiers, certificate and witness data, the specs they
  implement, and the paper's own prose mathematics.
- T2, supporting: outputs appear publicly but are re-derived or
  audited by a T1 artifact (generators, corroborating heuristics,
  sampled diagnostics explicitly labeled non-certifying).
- T3, exploratory: nothing public rests on it.

## Required evidence

T1: a written contract, even one page (the claim, the error model,
the environment assumptions, the governing soundness rule in
incapable-of-falsely-reporting form); adversarial review to closure,
including at least one independent, decorrelated review and a closure
verdict ON THE FINAL SOURCE as shipped; negative or broken-variant
tests wherever the artifact has refusal paths; runs recorded in at
least two environments, three where feasible; runtime
self-verification of every decimal/binary threshold; version banner;
in-file changelog.

T2: code review (any family); outputs cross-checked against a T1
artifact or an independent reimplementation; banner and changelog.

T3: changelog, and an explicit exploratory marking.

## Decorrelation

No T1 artifact closes on a single reviewer's evidence alone. An
independent, decorrelated review is mandatory; a second is required when
the artifact lacks a negative-test suite (the suite is the stronger
substitute).

## Machine-checked claims

"Machine-computed constants only" extends to reasoning ABOUT floats:
any claim of rounding direction, ulp distance, or threshold ordering
(decimal vs binary) is established by an executed exact-arithmetic
check, never by prose. Where the claim guards a verdict, the check
runs at runtime in the artifact itself.

## Review protocol

- Fresh-context review: the artifact, its data, and minimal
  ground truth only; in-file claims are audited as claims; prior
  findings are never disclosed (the artifact's own changelog is the
  only permitted disclosure).
- Closure round: final source as shipped, adjudications disclosed,
  explicit CLOSED/REOPEN verdict. A review of version N does not
  close version N+1.
- Stopping rule declared at dispatch, before results arrive. Default:
  a clean round (no CRITICAL or HIGH soundness finding) closes; an
  adopted soundness finding reopens -- fix, version bump, rerun,
  and only the finder re-reviews.
- Every finding is adjudicated individually and in writing: ADOPTED,
  NO DEFECT (with the demonstration), or DECLINED (with the reason).
  No-defect and declined findings do not block closure but are
  recorded.
- In-context self-checks are same-context checks, labeled as such;
  they never substitute for a fresh-context round.

## Mechanics (binding)

- Version banners on every executable artifact; a missing banner
  means a stale file.
- In-file changelogs; never silent edits. Ledgers are corrected by
  marked in-place correction plus a changelog entry, never rewritten.
- Discrepancies are adjudicated, never tolerated; the adjudication is
  recorded where the discrepancy lives.
- A verification source that fails ground truth is terminated, never
  laundered -- soundness defects are fixed regardless of margin.
- Frozen-in-flight: after gate closure an artifact changes only with
  a version bump, a rerun, and a re-recorded result.
- Environment assumptions (rounding model, libm error bounds,
  arithmetic semantics) are stated in the contract and recorded with
  the closure verdict.


## Relation to compute policy

`COMPUTE-POLICY-v1.md` governs performance rewrites, parallel scheduling,
runtime metadata, and benchmark reporting. For T1 artifacts, any optimization
that changes evaluator logic, arithmetic semantics, thresholds, split rules,
or certificate predicates reopens evidence obligations under this standard.
Pure scheduler changes still require sequential-vs-parallel validation. For
this repo COMPUTE-POLICY is otherwise N/A: every check is a small exact-`sympy`
or finite computation with no long-running, parallel, or budgeted compute, so
no binding `COMPUTE-POLICY-v1.md` copy is shipped.

## Instantiation for this repository

Artifact -> tier (`verification/STATUS-LEDGER.md` is the claim table; this maps
artifacts to evidence load):

- **T1 (certifying).** `code/exact_checks.py` (C1 Ptolemy, C3 sign-uniqueness
  exact, C4 pentagon/far-commute); `code/diagonal_form.py` (C8 exact identity
  `fᵀQf=Q` + composition); `code/exact_kernel.py` (C9 exact, inner `PB_{3,4}`);
  `code/appendix_full_matrix.py` (C10 full-matrix audit / erratum); the paper's
  prose mathematics (Lemma 1, Propositions 1–5, the Observation); and the
  independent `BLIND_CHECK.py`.
- **T2 (supporting).** `code/relations.py`, `code/sign_uniqueness.py`,
  `code/pentagon_typo.py`, `code/appendix_reconstruction.py`, `code/float_kernel.py`
  — floating-point smoke tests / reconstructions, each superseded and
  cross-checked by its T1 exact counterpart.
- **T3 (exploratory).** `code/kernel_probe.py` — the diagnostic
  braid → Delaunay-flip pipeline; used only to detect the combinatorial flip
  sequence, which `code/exact_kernel.py` then re-applies in exact arithmetic.

Decorrelation status: authored and then cross-checked by an independent,
decorrelated review and an independent hand re-derivation; the executed
suite has run in two environments. C2 and C10 previously carried a transcription
dependency on the source paper; it was discharged (2026-06-27) by a direct read
and an independent fresh-context re-read of arXiv:2606.20473v1, both confirming
the pentagon typo and the `I+skew` result matrices.

Changelog. 2026-06-12 v1 template copy created with the lab standard.
2026-06-26 binding copy added to this repo; Instantiation filled (artifact->tier).

#!/usr/bin/env bash
# Reproduces every claim in the note, in claim order. The exact certificate runs first;
# total runtime is a few minutes (exact_kernel.py dominates). No network access, no seeds needed.
set -e
echo "=== cbg-2026-representation-checks reproduction suite. Exact certificates run first; a few minutes total. Claim statuses and evidence types: verification/STATUS-LEDGER.md. ==="
# Resolve a working interpreter: python3 (Linux/macOS), else python (Windows).
PY=python3
"$PY" -c "" 2>/dev/null || PY=python
"$PY" -c "" 2>/dev/null || { echo "no working Python interpreter found"; exit 1; }
echo "== C1/C3/C4: exact certificate (Ptolemy, pentagon 120/120, sign-uniqueness, far-commute) =="
"$PY" code/exact_checks.py
echo "== C8: diagonal similarity f'=D f D^-1 (f^T Q f = Q; proof + confirmation) =="
"$PY" code/diagonal_form.py
echo "== C2: main-text pentagon typo =="
"$PY" code/pentagon_typo.py
echo "== C5/C6: convention lock + rational pure-braid property on b_4,7 =="
"$PY" code/appendix_reconstruction.py
echo "== C10: appendix worked-example erratum (full-matrix audit) =="
"$PY" code/appendix_full_matrix.py
echo "== C4: groupoid relations (fast smoke test) =="
"$PY" code/relations.py
echo "== C3: sign-uniqueness search + negative tests (fast smoke test) =="
"$PY" code/sign_uniqueness.py
echo "== Braid -> Delaunay-flip pipeline (diagnostic): braid relation, far-commute, Borromean -> I =="
"$PY" code/kernel_probe.py
echo "== C9: invariant = abelianization, ker f' = [PB,PB] (inner n=3,4; floating point) =="
"$PY" code/float_kernel.py
echo "== C9: exact lift (sympy; inner PB_3, PB_4) =="
"$PY" code/exact_kernel.py
echo "== C9 obligation (R2): exact/combinatorial flip-sequence certificate (Checker B: B1 legal-loop + B2-lite robustness) =="
"$PY" code/check_flip_sequence.py
echo "== Done: compare printed values against the README claim map. =="

# Author's independent check (before reading code/ or verification/)

Recompute the headline quantities from the operator formulas of the source paper (arXiv:2606.20473)
and this note, WITHOUT opening `code/`. Inputs: only `zeta_l = l` and the formulas listed here.

1. **The two printed appendix flip matrices** `A_4615`, `A_1524`. With the positional convention
   `(i,j,k,l) = (a,c,b,d)` and `cos = sqrt(z_il z_jk / (z_ik z_jl))`,
   `sin = sqrt(z_ij z_kl / (z_ik z_jl))` (`z_ab = zeta_a - zeta_b`), confirm
   `A_4615 -> cos 0.79056942, sin 0.61237244` and `A_1524 -> cos 1.06066017, sin 0.35355339 i`,
   matching the paper's printed values.
2. **Ptolemy/Plucker** `z_il z_jk + z_ij z_kl = z_ik z_jl`, hence `cos^2 + sin^2 = 1`, over every
   ordering.
3. **The single rational flip preserves the diagonal form** `Q(D_pqr) = 1/(z_pq z_pr z_qr)`:
   `f^T Q_p f = Q_t` (symbolic, exact).
4. **The worked pure braid** `b_4,7` (rational `f`, 16 named flips, rightmost acts first):
   determinant 1, unit diagonal, rank-2 one-step unipotent `((M-I)^2 = 0)`.
5. Only now open `code/` and `verification/` and compare. On disagreement: **stop and report; do not
   reconcile by editing the check to match.**

`BLIND_CHECK.py` implements (1)-(4) from scratch, importing nothing from `code/` (only `numpy`,
`sympy`, and the formulas above). It is deliberately a second, independent implementation of the
touched-triangle composition for step (4), so a transcription or convention error in `code/` would
surface as a disagreement here.

## Completion record

Run `python BLIND_CHECK.py` from the repository root. Expected: exit 0 and `BLIND CHECK PASSED`, with
the four blocks reproducing the printed flip matrices, `cos^2+sin^2=1` to ~1e-16, the exact
`f^T Q_p f - Q_t = 0`, and the `b_4,7` active block at `det 1`, `max|diag-1| ~ 1e-16`, `rank(M-I)=2`,
`(M-I)^2=0`. Record the date and outcome of your run here.

- Run date: _______  Outcome: PASS / FAIL (if FAIL, paste the failing block and stop).

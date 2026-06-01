# Pivot plan: the U_full row-sum assembly (PLANNED, NOT ACTIONED)

Written 2026-06-01 (Session 39) as a cold-start brief. Do not action yet; pick up when Opus
credits are available. **First action on resume: re-read this file, then `git fetch` and rebase the
working branch on the latest base before doing anything else.**

## Where we are (what is already proved)

- **Lemma A (clean), all k:** `||P_a U_clean P_b||_2 = 2^{-(b-a)/2}` exact, for all `0 <= a < b <= k-2`.
  Proof + referee record in
  [HALFSHIFT_S4_LEMMA_A_PROOF.md](HALFSHIFT_S4_LEMMA_A_PROOF.md).
  Reduces to two elementary 2-group lemmas (CU, SB) + the trivial geometric series S4.
- **Half-Shift Invariance / coset-uniformity:** proved (same doc, Sections 1-2), so the foundation's
  R1 (`U_clean` strict-upper) and R2 (rank-1 defect `D = e_{r*} c^*`) are unconditional.
- **Lemma B (collision bound), all k:** `||tril(Q_k)||_2 <= sqrt(3) 2^{-k/2} -> 0`, unconditional.
  See [LEMMA_B_PROOF.md](LEMMA_B_PROOF.md).

The certificate runs on `U_full = U_clean + D`. `Q_k[a,b] := ||P_a U_full P_b||_2`.

## The one remaining math obligation

> Prove, for all `k`: `max_a  sum_b  Q_k[a,b] * 2^{a-b}  <  1`  (equivalently `rho(Q_k) < 1`),
> hence `|lambda_2(U_k)| < 1` uniformly, hence **no non-trivial Syracuse cycles**.

This is the certificate-level assembly of Lemma A (upper) + Lemma B (defect). It is NOT a trivial
wrap-up; see the crux below.

## What is in hand for the assembly

- **Upper part `a < b`:** `Q[a,b] ~ 2^{-(b-a)/2}` (Lemma A, exact for `U_clean`; `U_full` adds an
  `O(2^{-k})` leak, measured). Its weighted row-sum contribution is the exact geometric series
  `sum_{d>=1} 2^{-d/2} 2^{-d} = 2^{-3/2}/(1-2^{-3/2}) = 0.5469`, constant in `k`.
- **Lower+diagonal part `a >= b`:** entirely the defect (R1). `Q_D[a,b] = u_a v_b` rank-1, with
  `u_a = 2^{-(a+1)/2}` (since `u_a^2 = d_a/N = 2^{-1-a}`) and `v_b = ||P_b c||`, `sum_b v_b^2 = ||c||^2
  <= 3 * 2^{-k}`.
- **Measured ground truth (build_T):** row-sum flat `~0.6345` for `k=6..12`; `rho(Q) ~ 0.566`; upper
  part `0.5469` exact; lower+diag part `~0.156` and decaying. Margin to 1 is `0.366`.

## THE CRUX of the assembly (do not hand-wave this)

The lower+diagonal weighted row-sum for row `a` is
```
  S_a := sum_{b <= a} Q_D[a,b] 2^{a-b} = u_a 2^a sum_{b<=a} v_b 2^{-b}
       = 2^{(a-1)/2} * sum_{b<=a} v_b 2^{-b}.
```
The prefactor `2^{(a-1)/2}` GROWS with `a` (up to `~2^{k/2}` at the top row `a = k-2`). For `S_a` to
stay bounded (and `< 0.453` so the total is `< 1`) uniformly in `k`, the `v_b` profile must decay fast
enough to beat the `2^{(a-1)/2}` weight. The measured `S_a` is flat `~0.156`, so it DOES hold -- but
the analytic mechanism is the real content and has not been written.

**The lever:** the `v_b` profile. From the defect covector `c` (the partial Gauss sum at `r*`), the
numerics (`adv_tril_sep_correct.py`, `_reconcile_out.txt`) show the "gamma profile"
`gamma_b * 2^b * 2^{k/2} -> ~1.56`, i.e. roughly `v_b ~ C * 2^{b - k/2}`. Plugging in:
`sum_{b<=a} v_b 2^{-b} ~ C 2^{-k/2} sum_{b<=a} 1 ~ C a 2^{-k/2}`, so
`S_a ~ 2^{(a-1)/2} * C a 2^{-k/2}`. At `a ~ k`: `~ 2^{k/2} * C k * 2^{-k/2} = C k` -- which GROWS.
The flat measured value means this crude estimate is wrong: the `v_b` must decay (not be flat in `b`)
near the top, OR the sum is dominated by small `b` where `v_b 2^{-b}` is largest. **Pin the exact `v_b`
profile** (it has a closed form via the `r*` defect / the `gamma` data) and bound `S_a` correctly.
This weight-vs-decay battle is the genuine remaining work -- treat it as a third lemma, not a wrap-up.

- Also handle the **upper-leak**: `U_full`'s upper entries = `U_clean` (exact `2^{-d/2}`) + `O(2^{-k})`
  defect leak. Show the leak's weighted-row-sum contribution stays summable and `-> 0`, so the upper
  bound stays `0.5469 + o(1)`.

## Suggested method (the writer + adversarial-referee pattern that has worked all session)

1. **Pin the `v_b` profile** in closed form from the `r*` defect covector (start from `adv_tril_sep_correct.py`'s
   `level_energies` and the `gamma` profile in `_reconcile_out.txt`).
2. **Writer:** prove `max_a S_a < 0.453` (ideally `-> 0`) for all `k`, plus the upper-leak bound, then
   assemble `row-sum(Q_full) <= 0.5469 + o(1) < 1`.
3. **Adversarial referee (workflow):** attack the bound at the row `a` that maximises `S_a`; both
   parities; check the `2^{a-b}` weight vs `v_b` decay at the top row `a=k-2`; verify the analytic
   bound dominates the measured `0.156` and stays `< 0.453` for `k` pushed well past 12; hunt for a
   Lemma-B-style boundary/parity surprise.
4. If it holds: `|lambda_2(U_k)| < 1` uniformly => no non-trivial Syracuse cycles. Then **Lean
   formalisation** (CU, SB, S4, the counting, Lemma B's collision bound, the assembly) and a paper.

## Oracles / files for the next session

- `HALFSHIFT_S4_LEMMA_A_PROOF.md`, `LEMMA_B_PROOF.md` -- the two proved lemmas.
- `adv_tril_sep_correct.py` -- `v_b = ||P_b c||`, `tril(Q_D)`, the defect covector; the `v_b` profile lives here.
- `audit_halfshift_s4.py`, `attack1_lemmaA_proof.py` -- Lemma A / CU / S4 verification.
- `_reconcile_out.txt` -- the gamma profile and the measured certificate (`row-sum 0.634`, `rho 0.566`).
- `analytic_proofs.py` : `build_T(k)` -- the ground-truth `Q_full` row-sum oracle.
- Targets to reproduce: upper `0.5469` exact, lower+diag `~0.156 -> 0`, total `~0.6345`, `rho ~0.566`.

## Definition of done

All-`k` proof that `rho(Q_k) < 1` (the assembly above). Then the certificate proves **cycle
elimination** unconditionally. Then Lean + paper.

## Honest ceiling (keep this in every writeup)

Cycle elimination is **not** Collatz. The divergent-trajectory half is untouched and is the
irreducibility statement (Mori 2024) equivalent to the full conjecture. This program, complete, would
be a hard, publishable "no non-trivial Collatz cycles," strictly weaker than the conjecture.

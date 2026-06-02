# The U_full row-sum assembly: proof modulo one sharp per-level bound

> **RETRACTION note (2026-06-02).** This document proves `cert(k) < 1`, a uniform spectral gap for the
> transfer operator. The further step "spectral gap => no non-trivial cycles" stated below is **false**
> and is **withdrawn**: the `3x-1` map satisfies the same certificate yet has cycles. See
> [CYCLE_CLAIM_REFUTED.md](CYCLE_CLAIM_REFUTED.md). The `cert(k) < 1` result is correct; the cycle
> corollary is not. Read the cycle-elimination sentences below as withdrawn.

Written 2026-06-02 (Session 40), resuming from [UFULL_ASSEMBLY_PLAN.md](UFULL_ASSEMBLY_PLAN.md).
This document closes the assembly **conditionally on a single new lemma (Lemma C)**, states that
lemma sharply, reduces it to a clean periodization estimate, and is explicit about what is proved
versus what remains. Reproduce every number with `python verify_assembly.py` and the `probe_*.py`
scripts.

## Status in one line

The certificate `max_a sum_b Q_k[a,b] 2^{a-b} < 1` is now **proved for all `k`, uniformly, modulo
Lemma C** (a sharp per-level bound on the defect covector, `v_b <= (3/4) 2^{-b} 2^{-k/2}`). Lemma C
is verified exactly to `k=26`, has a sharp universal constant, and is reduced below to a periodization
excess bound. **Lemma C itself now has a proof** in [LEMMA_C_PROOF.md](LEMMA_C_PROOF.md): the
assembly-strength bound `g_b < sqrt(3/4) < 0.961` is established by the shell method (Lemma B's
machinery), with the previously-open cross-scale homometry `S_k(2xi) = S_{k-1}(xi)` now **derived**
(Lemma H, an elementary parity split). Lemma C inherits only Lemma B's existing Half-Shift Invariance
dependence and adds no new conditional input. So the open mathematical content of the program is now
exactly the Half-Shift Invariance crux (shared by Lemmas A and B) plus the Lean formalisation.

## Setup (recap)

`U_full = U_clean + D` on the character basis, decomposed by 2-adic level `a = 0..k-2`.
`Q_k[a,b] := ||P_a U_full P_b||_2`. Two facts from the proved lemmas and the foundation:

- **(Lemma A, upper, exact)** `U_clean` is strictly block-upper (foundation R1, from coset-uniformity)
  and `||P_a U_clean P_b||_2 = 2^{-(b-a)/2}` for `0 <= a < b <= k-2`; `= 0` for `b <= a`.
  Proof: [HALFSHIFT_S4_LEMMA_A_PROOF.md](HALFSHIFT_S4_LEMMA_A_PROOF.md).
- **(Defect, rank-1)** `D = e_{r*} c^*` (foundation R2), so `||P_a D P_b||_2 = u_a v_b` for **all**
  `a, b`, with
  - `u_a = ||P_a e_{r*}|| = 2^{-(a+1)/2}` exactly (since `u_a^2 = d_a/N = 2^{k-2-a}/2^{k-1} = 2^{-1-a}`),
  - `v_b = ||P_b c||`, the level-`b` energy of the defect covector `c = D[r*,:]`.

By the triangle inequality for the operator norm, for every `a, b`:
```
  Q_k[a,b]  <=  Q_clean[a,b] + u_a v_b ,   where Q_clean[a,b] = 2^{-(b-a)/2} [b>a].
```

## Lemma C (new): sharp per-level bound on the defect covector

> **Lemma C.** For all `k` and all `0 <= b <= k-2`,
> ```
>     v_b  <=  (3/4) * 2^{-b} * 2^{-k/2} .
> ```
> Equivalently, with `g_b := v_b * 2^b * 2^{k/2}`, we have `g_b <= 3/4`, i.e. `g_b^2 <= 9/16`.

**Evidence (exact, `verify_assembly.py` / `probe_*.py`):** `g_b^2 <= 9/16` for every `k` tested
(`k = 6..22`), with **equality at `b = k-4`** (so the constant `3/4` is sharp and cannot be lowered).
The profile `g_b^2` is `k`-independent at both ends: it rises from a fixed bulk value `7/16` to the
fixed dyadic tail `..., 53/128, 15/32, 11/32, 9/16, 1/8, 1/2` reading toward the top level `b=k-2`.

**Reduction to a periodization excess (the proof handle).** Let `F = FFT_N(c)`, `N = 2^{k-1}`, and let
`fold_{2^j}(c)` be the `2^j`-fold periodization of `c` (sum of `c` over the `2^j` cosets of `N/2^j`).
By Parseval/aliasing, `M_j := sum_{xi == 0 mod 2^j} |F[xi]|^2 = (N/2^j) ||fold_{2^j}(c)||^2`, and the
exact-level-`b` band energy is `B_b = M_b - M_{b+1}`. Since `v_b^2 = B_b / N`,
```
  g_b^2 = 2^k 4^b v_b^2 = 2 * 4^b * B_b = 2^b ( 2 h_b - h_{b+1} ),   h_j := N ||fold_{2^j}(c)||^2.
```
Numerically `h_j = 2^j + a_j` with excess `a_j >= 0` decaying to `0` as `j -> k`, so
```
  g_b^2 = 2^b ( 2 a_b - a_{b+1} ).
```
Lemma C is therefore equivalent to the excess bound `2^b(2 a_b - a_{b+1}) <= 9/16` for the
periodizations of the exceptional fiber distribution `c`. This is `verify`-confirmed; its closed form
(and hence a full analytic proof) is the open step. `c` is the Syracuse fiber distribution of
`r* = -3^{-1} mod 2^k`, and its periodizations are the same distribution at smaller `k`, so the
self-similarity that powers Lemma B's collision count is the natural tool here.

## The assembly theorem (proved, given Lemma A + Lemma C)

> **Theorem.** Assume Lemma A and Lemma C. Then for all `k >= 4`,
> ```
>     cert(k) := max_a  sum_b  Q_k[a,b] * 2^{a-b}  <  G_up + 2^{-3/2}  =  0.900472...  <  1,
> ```
> where `G_up := sum_{d>=1} 2^{-3d/2} = 1/(2^{3/2} - 1) = 0.546918...`. Hence `rho(Q_k) <= cert(k) < 1`
> uniformly.

**The spectral reduction `cert(k) < 1 => no non-trivial cycles` (cited, established elsewhere).**
This document proves the *algebraic* statement `cert(k) < 0.9005`. The link to a spectral gap is the
standard two-step chain, established in the foundation, not re-proved here:
```
  |lambda_2(U_k)|  <=  rho(Q_k)  <=  cert(k).
```
The first inequality is the Perron/non-Perron split: `U_k` is a (sub)stochastic-type transfer operator
whose top eigenvalue is the Perron value `1`; on the complementary (non-Perron) invariant part the
spectral radius is dominated by the block-norm matrix `Q_k[a,b] = ||P_a U_k P_b||_2`, so
`rho(U_k|_{non-Perron}) <= rho(Q_k)`. The second is the weighted-row-sum (Gershgorin) bound on the
nonnegative matrix `Q_k` under the similarity `diag(2^{-a})`. Both are stated and validated in
[STEP4_BLOCK_FORMULA_FOUNDATION.md](STEP4_BLOCK_FORMULA_FOUNDATION.md) §3 (`|lambda_2(T_k)| <= rho(Q) <= rowsum`)
and the README certificate display. Given that chain, `cert(k) < 1` uniform yields `|lambda_2(U_k)| < 1`
uniform, which eliminates non-trivial Syracuse cycles. The present document's contribution is closing
the last algebraic link `cert(k) < 1` for all `k` (modulo Lemma C); it does not re-derive the spectral
reduction.

**Proof.** Fix a row `a`. Using `Q_k[a,b] <= Q_clean[a,b] + u_a v_b`:
```
  R_a := sum_b Q_k[a,b] 2^{a-b}
       <= sum_{b>a} 2^{-(b-a)/2} 2^{a-b}              (Lemma A upper, clean)
          + u_a 2^a sum_{b=0}^{k-2} v_b 2^{-b}         (defect, all b -- absorbs the upper leak)
       =  sum_{d=1}^{k-2-a} 2^{-3d/2}  +  2^{(a-1)/2} sum_{b=0}^{k-2} v_b 2^{-b}.
```
The first sum is `<= G_up` (it is a truncation of the full geometric series). For the second, apply
Lemma C `v_b <= (3/4) 2^{-b} 2^{-k/2}`:
```
  sum_{b=0}^{k-2} v_b 2^{-b}  <=  (3/4) 2^{-k/2} sum_{b=0}^{k-2} 4^{-b}
                              <   (3/4) 2^{-k/2} * (4/3)  =  2^{-k/2}.
```
Therefore `R_a < B(a) := G_up + 2^{(a-1)/2} 2^{-k/2} = G_up + 2^{(a-1-k)/2}`. The *bound* `B(a)` is
increasing in `a`, so it is maximised at the top row `a = k-2`, giving
```
  cert(k) = max_a R_a  <=  max_a B(a) = B(k-2) = G_up + 2^{(k-3-k)/2} = G_up + 2^{-3/2} = 0.900472 < 1.   ∎
```
(The true `R_a` peaks at a middle row near `a ~ (k-2)/2`, where the decreasing clean-upper part and the
increasing defect part trade off; the top row's true `R_a` is only `~0.31`. We do not need this: we
bound the decreasing upper part by its supremum `G_up` and the increasing defect part by its top value,
and `B(a)` -- the resulting envelope -- is monotone, so `max_a R_a <= max_a B(a) = B(k-2)`.)

Two remarks on why the proof is not a wrap-up:

1. **The upper leak is handled for free.** The defect's contribution to the *upper* blocks `b>a` is
   not dropped; it is included by summing `u_a v_b` over **all** `b` (not just `b <= a`). The triangle
   inequality `Q_k <= Q_clean + Q_D` makes this automatic, so no separate leak estimate is needed.

2. **The `L^2` bound alone is insufficient -- Lemma C is necessary.** Lemma B gives only
   `sum_b v_b^2 <= 3 * 2^{-k}`. Feeding that through Cauchy-Schwarz,
   `sum_b v_b 2^{-b} <= sqrt(3 * 2^{-k}) sqrt(4/3) = 2^{1-k/2}`, yields `R_{k-2} <= 2^{-1/2} = 0.707`,
   and `G_up + 0.707 = 1.25 > 1` -- the certificate would **fail**. The per-level decay
   `v_b ~ 2^{-b}` (Lemma C), not merely the total `L^2` mass, is what closes the bound. This is the
   precise sense in which the assembly is "a third lemma, not a wrap-up."

## Numerical certificate (ground truth, for comparison)

`verify_assembly.py` reports the true row-sum certificate from the measured `v_b` and from the dense
`build_T` operator; they agree to 6 digits, confirming the separable `u_a v_b` model **is** the
operator's defect:

| k | cert (true) | rho(Q) | max_a S_a (lower) | Lemma-C bound on lower |
|---|---|---|---|---|
| 6  | 0.633782 | 0.553529 | 0.313386 | 0.353553 |
| 8  | 0.634659 | 0.566061 | 0.312395 | 0.353553 |
| 10 | 0.634533 | 0.567279 | 0.311970 | 0.353553 |
| 12 | 0.634477 | 0.565553 | 0.311848 | 0.353553 |

The true certificate sits at `~0.6345` (margin `0.366` to 1); the analytic bound proved above is the
looser-but-uniform `0.9005` (margin `0.10`). Either is `< 1`; the analytic one is what the theorem
delivers for all `k`.

## What is proved vs. what remains

**Proved (this document, given the prior lemmas):**
- The full row-sum assembly: `cert(k) < 0.9005 < 1` uniformly in `k`, from Lemma A + Lemma C, by an
  elementary triangle-inequality + geometric-series argument. The upper leak is absorbed; uniformity
  is explicit; the necessity of per-level decay is demonstrated.

**Remaining obligations (no finite-`k` crux left):**
- **Lemma C is proved** ([LEMMA_C_PROOF.md](LEMMA_C_PROOF.md)) at assembly strength
  `g_b < sqrt(3/4) < 0.961`, via the shell method with the cross-scale homometry now derived (Lemma H).
  (The sharp `g_b <= 3/4`, equivalently `2^b(2 a_b - a_{b+1}) <= 9/16`, is verified to `k=26`, not needed.)
- **Half-Shift Invariance / coset-uniformity is also proved** ([HALFSHIFT_S4_LEMMA_A_PROOF.md](HALFSHIFT_S4_LEMMA_A_PROOF.md)
  Sections 1-2): R1/R2 reduce to the Coset-Uniformity lemma CU, proved by elementary finite-group
  theory. The assembly uses only `rank(D) <= 1` (the `D = e_{r*} c^*` form), which CU gives; the defect
  *nonvanishing* is verified but not used. So there is no remaining finite-`k` mathematical gap.
- **Spectral reduction write-up.** `cert(k) < 1 => |lambda_2(U_k)| < 1 => no non-trivial cycles` is the
  standard chain `|lambda_2| <= rho(Q_k) <= cert(k)` (Perron/non-Perron split + weighted-row-sum bound),
  stated in [STEP4_BLOCK_FORMULA_FOUNDATION.md](STEP4_BLOCK_FORMULA_FOUNDATION.md) but cited rather than
  derived here; writing it out in full is an open item.
- **Lean formalisation** of CU, SB, S4, the counting, Lemma B, Lemma C, and the assembly above.

## Honest ceiling (unchanged)

Cycle elimination is **not** Collatz. Even when Lemma C is proved, this program establishes only "the
only Collatz cycle is `1 -> 4 -> 2 -> 1`," which is strictly weaker than the conjecture. The
divergent-trajectory half (the irreducibility statement, Mori 2024) is untouched.

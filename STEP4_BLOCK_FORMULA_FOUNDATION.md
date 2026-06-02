# Step-4 Block-Formula Foundation (2026-06-01)

> **RETRACTION note (2026-06-02).** The block formula and `Q[a,b]` structure here are correct and feed
> the operator spectral-gap proof. The implied chain `|lambda_2(T_k)| < 1 => no non-trivial cycles` is
> **false and withdrawn**: the `3x-1` map gives the same `|lambda_2| <= rho(Q) < 1` and has cycles. See
> [CYCLE_CLAIM_REFUTED.md](CYCLE_CLAIM_REFUTED.md). The spectral-gap content stands; the cycle corollary
> does not.

Exact, validated foundation for the step-4 entry-bound lemma of the uniform spectral gap.
Everything here is checked against `build_T` to machine precision; finite-k verification is named
as such and never conflated with all-k proof.

Scripts (all under [Collatz](Collatz)):
- [diagnose_block_formula.py](diagnose_block_formula.py) - corrects the block formula, proves the Gauss-sum collapse.
- [step4_foundation.py](step4_foundation.py) - U = U_clean + rank-1 defect; Q tables.
- [step4_within_level.py](step4_within_level.py) - the isometry identity B*B = 2^-d I.
- [step4_proof_mechanism.py](step4_proof_mechanism.py) - within-level mechanism + defect covector.
- [step4_deliverable.py](step4_deliverable.py) - consolidated deliverable tables.

## 0. CORRECTION to the stated block formula (this matters - it was a trap)

The prompt's literal formula
```
(U chi_eta)(r) = [v(r) <= b] w^{eta(3r+1)/2^{v(r)}},   b = v2(eta), v(r)=v2(3r+1)
```
does NOT match `build_T` (`max|U chi - formula| ~ 0.25`, flat in k). The honest exact identity,
verified to 3e-16, has TWO additions:

```
(U chi_eta)(r) = [v(r) < k] [v(r) <= b]  w^{ eta (3r+1)/2^{v(r)} }        for v(r) < k
(U chi_eta)(r*) = (1/2^k) sum_{m<2^k} w^{ eta Syr(r* + m 2^k) }            for v(r*) >= k
```

- The `v(r) < k` guard is essential: the coset-collapse (Gauss sum) is valid only when the Syracuse
  image of `r + m 2^k` sweeps a FULL coset of size `2^v` with uniform multiplicity `2^{k-v}`. That
  holds iff `v(r) < k` (so `k - v >= 1`, the half-shift/coset-uniformity lemma applies).
- Exactly ONE residue fails it: `r* = (4^{ceil(k/2)}-1)/3`, the unique odd `r` with `v2(3r+1) >= k`
  (the same `r*` as the rank-1 S_odd theorem). At `r*` the value is a genuine partial Gauss sum.

So  **U = U_clean + D**  where U_clean is the masked-phase operator (rows `v(r) < k`) and
**D is rank-1, supported on the single row `r*`**.

### The Gauss-sum collapse (the heart of U_clean), derived honestly
For `v(r) < k`, coset-uniformity gives
```
(U chi_eta)(r) = (1/2^v) sum_{j=0}^{2^v-1} w^{ eta ( (3r+1)/2^v + j 2^{k-v} ) }
              = w^{eta(3r+1)/2^v} * (1/2^v) sum_{j} W^j ,   W = w^{eta 2^{k-v}} = exp(2 pi i eta / 2^v).
```
The inner geometric sum is `2^v` if `2^v | eta` (i.e. `v2(eta) = b >= v`) and `0` otherwise
(`W^{2^v} = exp(2 pi i eta) = 1` always, so it is a clean indicator). Hence
```
(U chi_eta)(r) = [v <= b] w^{eta(3r+1)/2^v}     (v = v(r) < k).
```
This is the masking `[v(r) <= b]` - it is the Gauss-sum annihilation of source frequencies coarser
than the local valuation. PROVED for all k given coset-uniformity (Half-Shift Invariance Lemma).

## 1. Q[a,b] = ||P_a U P_b||_2 in closed form

Decompose `g_eta(r) = [v(r)<k][v(r)<=b] w^{eta q(r)}`, `q(r) = (3r+1)/2^{v(r)}`, over levels.
The map `eta -> q(r)`-phase preserves frequency level (x3 is a 2-adic unit; +1 is the carry). The
block `B = P_a U_clean P_b` (target level a, source level b) has the EXACT structure:

**Upper (a < b, d = b - a >= 1):  B is `2^{-d/2}` times an ISOMETRY of level b into level a.**
Verified `B* B = 2^{-d} I_{d_b}` to ~1e-10 (k<=12), where `d_b = 2^{k-2-b}` is the source dim.
- All `d_b` singular values equal `2^{-d/2}` exactly (spread ~1e-10).
- => blocks are **NOT rank-1**; they are full-source-rank scaled isometries.
- => `||B||_2 = 2^{-d/2}` EXACT (this is Lemma A's value, for U_clean exactly; for U_full to
  `O(2^{-k/2})`, since the only difference is the rank-1 defect row).

**Diagonal/lower (a >= b):  P_a U_clean P_b = 0  EXACTLY** (verified to 5e-11).
U_clean is strictly upper-triangular in level. So ALL diagonal+lower entries of Q_full come purely
from the rank-1 defect D. `||tril(Q_full)||_2 = beta_k 2^{-k/2}`, `beta_k -> ~0.684`.

Numerical confirmation (`step4_deliverable.py`, validated vs build_T):
```
k | max|Q_up - 2^(-d/2)| | ||tril Q||_2 (x 2^(k/2)) | rho(Q) | rowsum(D=2^j)
6 |      4.7e-04         |   0.0863  (0.6902)       | 0.5535 |   0.6338
8 |      1.4e-04         |   0.0429  (0.6860)       | 0.5661 |   0.6347
10|      3.7e-05         |   0.0214  (0.6849)       | 0.5673 |   0.6345
12|      9.4e-06         |   0.0107  (0.6846)       | 0.5656 |   0.6345
```
The upper deviation is `O(2^{-k/2})` (halves per k) - exactly the defect leaking into the upper
entries of Q_full. For U_clean the upper entries are `2^{-d/2}` to machine precision.

### Energy <-> norm reconciliation (Q_mass = 2^-d  vs  Q_norm = 2^-d/2)
`Q_mass[a,b] = (1/d_b)||B||_F^2`. For the isometry block: `||B||_F^2 = sum_{i=1}^{d_b}(2^{-d/2})^2
= d_b 2^{-d}`, so `Q_mass = 2^{-d}` EXACT. `Q_norm = 2^{-d/2}` is the single worst output direction.
No contradiction: the input energy spreads over `d_b` orthonormal output directions each of amplitude
`2^{-d/2}`; the operator norm sees one direction (`2^{-d/2}`), the mass sees the average (`2^{-d}`).
`Q_norm = sqrt(Q_mass)` because the block is an isometry up to the scalar `2^{-d/2}` - it is NOT the
`sqrt(d_b * Q_mass)` rank-1 bookkeeping suggested in the prompt's hint (that would only hold if the
block were rank-1, which it is not).

## 2. EXACT statements of the two lemmas that suffice for the certificate

**Lemma A (upper cascade).** For all k and all `0 <= a < b <= k-2`,
```
|| P_a U_k P_b ||_2  <=  2^{-(b-a)/2}.
```
What it needs: (i) the Gauss-sum collapse identity `(U chi_eta)(r) = [v(r)<=b] w^{eta q(r)}` for
`v(r) < k` (PROVED via coset-uniformity / Half-Shift Invariance), and (ii) the WITHIN-LEVEL isometry
identity `B*B = 2^{-(b-a)} I_{d_b}`. Identity (ii) is the one structural fact still to be proved for
all k. Its content: the level-a Fourier coefficients of `{g_eta : v2(eta)=b}` are orthonormal up to
the scalar `2^{-(b-a)/2}` - i.e. the carry map `r -> q(r)` reindexes source level-b characters into
level-a characters injectively with flat amplitude. (For U_full add the `O(2^{-k/2})` defect: Lemma A
for U_full holds with an additive `O(2^{-k/2})` slack, which the certificate absorbs.)

**Lemma B (diagonal + lower back-flow).** For all k,
```
|| (diagonal + lower part of Q_k) ||_2  <=  beta 2^{-k/2},   beta universal (measured ~0.684 -> 0.764 cap).
```
What it needs: (i) `P_a U_clean P_b = 0` for `a >= b` (PROVED: U_clean strictly upper in level, a
direct corollary of the masking `[v(r) <= b]` and `v(r) >= 1`), so the lower+diag block IS the
defect's level-matrix; (ii) the defect `D = e_{r*} c^*` is rank-1 with `||c||_2 = beta' 2^{-k/2}`
(`beta' -> 1.607`), and its level-matrix has norm `<= 0.764 2^{-k/2}`. The defect covector `c` is the
single partial Gauss sum at `r*`; its `2^{-k/2}` scale is the SAME scale as the rank-1 S_odd
inter-scale defect (Half-Shift Invariance), confirming the two phenomena coincide.
This is the cleaner alternative form (replaces the prompt's `beta 2^{-k/2} 2^{-(a-b)/2}` entrywise
claim, which is not needed - the rank-1 norm bound is enough and is exact).

## 3. Certificate closes (validated, finite k)

Given A and B, Step-3 row sum `max_a sum_b Q[a,b] 2^{a-b}`:
- upper part `<= sum_{d>=1} 2^{-d/2} 2^{-d} = 2^{-3/2}/(1-2^{-3/2}) = 0.547` (constant, all k);
- lower part `<= ||tril Q||_2 * (geometric scaling factor)`, bounded because `||tril Q|| = O(2^{-k/2})`
  and the diagonal scaling `2^{a-b}` over the at-most-(k-1) lower band contributes `<= O(2^{k/2})`,
  which the `2^{-k/2}` prefactor exactly cancels.
- Measured row-sum FLAT at 0.634 for k=6..13; `rho(Q) ~ 0.566`. Both `< 1` uniformly. Certificate
  closes => `|lambda_2(T_k)| <= rho(Q) <= 0.634`, conditional on A and B holding for all k.

## Honest status (PROVED-for-all-k vs VERIFIED-to-finite-k)

PROVED for all k (given the already-established Half-Shift / coset-uniformity lemma):
- The corrected exact block formula `U = U_clean + D`, masking `[v(r)<=b]` for `v(r)<k`.
- `U_clean` strictly upper-triangular in level (Lemma B part (i)).
- `D` rank-1, supported on the single row `r*` (Lemma B part (ii) structure).

VERIFIED to k=13, NOT yet proved for all k (the two remaining obligations):
- Lemma A identity (ii): `B*B = 2^{-(b-a)} I_{d_b}` (the within-level isometry / flat-amplitude fact).
  This is the crux. Route: show the carry map `r -> q(r)=(3r+1)/2^{v(r)}` induces, on level-b
  characters, an orthonormal family of level-a characters with constant amplitude `2^{-(b-a)/2}` -
  a Toeplitz/circulant secondary-character computation, plausibly a finite Gauss-sum identity.
- Lemma B scale: `||c||_2 = beta' 2^{-k/2}` with `beta'` bounded (measured `1.607`, flat). Route:
  bound the single partial Gauss sum at `r*` (connect to the S_odd defect column closed form).

## Obstruction note (a precise one beats a fake proof)

Lemma A is NOT provable by "the block is rank-1 so `||.||_2 = ||.||_F`" (the prompt's first hinted
route). The upper blocks have rank `d_b` (full source), NOT rank 1; their `d_b` singular values are
all equal. The correct lever is the **isometry identity `B*B = 2^{-d} I`**, a stronger statement than
rank-1, and the energy bookkeeping is `||B||_F^2 = d_b 2^{-d}` (NOT `||B||_2^2 = d_b Q_mass`). Any
attack that assumes rank-1 upper blocks will produce wrong dimension factors and must be rejected.

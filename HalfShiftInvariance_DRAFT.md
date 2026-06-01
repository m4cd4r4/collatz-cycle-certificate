# Half-Shift Invariance and the Rank-1 S_odd Theorem - Lean formalisation draft

Proven 2026-06-01 (Session 39). This replaces the hand-waved cancellation step in
analytic_proofs.py Lemma 1/2 with a correct, all-k proof. Target file: a new
`Collatz/HalfShiftInvariance.lean`, imported where the rank-1 fact is currently asserted.

## Definitions (Lean signatures)

```lean
-- Syracuse map on positive integers (already in Basic.lean as `syracuse` / `Syr`)
-- v2 = 2-adic valuation (already `padicValNat 2` or project `v2`)

-- Empirical transition probability at resolution k+1:
-- p x o = (1 / 2^(k+1)) * #{ m < 2^(k+1) | Syr (x + m * 2^(k+1)) ≡ o [MOD 2^(k+1)] }
-- (a rational; in Lean use a Finset.card / cast, or a counting measure)

-- Half-imbalance:
-- D x u = p x (2u+1) - p x (2u+1 + 2^k)
```

## Lemma 1 (Half-Shift Invariance) - the heart

```lean
theorem half_shift_invariance (k x : ℕ) (hx : Odd x)
    (hv : v2 (3*x+1) < k) :
    ∀ u, D k x u = 0
```

**Proof obligations (each is elementary):**

1. `factor`: write `3*x+1 = 2^v * q0` with `q0` odd, `v = v2(3*x+1)`.
2. `step_valuation`: for `n = x + m*2^(k+1)`, `3*n+1 = 2^v*(q0 + 3*m*2^(k+1-v))`, and since
   `k+1-v ≥ 2` the bracket is odd, so `v2(3*n+1) = v` and `Syr n = q0 + 3*m*2^(k+1-v)`.
   (Key arithmetic: `k+1-v ≥ 2` from `v < k`.)
3. `coset_uniform`: as `m` ranges over `ZMod (2^(k+1))`, `3*m` ranges over a complete residue
   system mod `2^v` (because `3` is a unit mod `2^v`, `Nat.Coprime 3 (2^v)`). Hence `Syr n mod 2^(k+1)`
   is uniform on the coset `q0 + (2^(k+1-v)) * ZMod (2^(k+1))`, with equal multiplicity `2^(k+1-v)`.
4. `shift_in_subgroup`: `2^k = 2^(k+1-v) * 2^(v-1)`, and `v ≥ 1` (since `3*x+1` is even for odd `x`),
   so `2^k` lies in the subgroup `2^(k+1-v) * ZMod (2^(k+1))`.
5. `translation_invariance`: a uniform measure on a coset is invariant under translation by a
   subgroup element, so `p x o = p x (o + 2^k)`, i.e. `D x u = 0`. ∎

Mathlib support: `ZMod`, `Nat.Coprime`, `ZMod.unitOfCoprime`, `AddSubgroup`, `Finset.card` counting,
`Nat.sub_add_cancel`. Step 3 is the only nontrivial one (a counting/uniformity argument); it can be
done by exhibiting the bijection `m ↦ 3*m` on `ZMod (2^v)`.

## Lemma 2 (uniqueness of the special residue)

```lean
theorem unique_special_residue (k : ℕ) (hk : 3 ≤ k) :
    ∃! r, r < 2^k ∧ Odd r ∧ 2^k ∣ (3*r+1)
-- and that r equals (4^(⌈k/2⌉) - 1)/3, with 3*r+1 = 2^(2*⌈k/2⌉)
```

**Proof:** `3` is a unit mod `2^k` (`Nat.Coprime 3 (2^k)`), so `3*r ≡ -1 [MOD 2^k]` has a unique
solution `r ≡ -3⁻¹`. It is automatically odd. The closed form `r* = (4^j-1)/3`, `j = ⌈k/2⌉`,
satisfies `3*r*+1 = 4^j = 2^(2j)` with `2j ≥ k`. (Already partially in RstarHierarchy.lean as
`rstar_formula`: `3 * rstar j + 1 = 4^j`.) ∎

## Theorem (Rank-1 S_odd)

```lean
theorem sodd_rank_one (k : ℕ) (hk : 3 ≤ k) :
    Matrix.rank (S_odd k) = 1
```

**Proof:**
- For every odd `r < 2^k` with `r ≠ r*`: `v2(3*r+1) < k`, and also `v2(3*(r+2^k)+1) < k`
  (adding `3*2^k`, which has valuation `k`, does not change bits below `k`). Lemma 1 gives
  `D_r = D_{r+2^k} = 0`, so the `S_odd` column at `r` is exactly zero.
- For `r = r*`: exactly one of the two lifts `{r*, r*+2^k}` mod `2^(k+1)` satisfies `2^(k+1) ∣ 3x+1`
  (Lemma 2 at level `k+1`), so it has `v2 ≥ k+1` and breaks half-shift invariance; the other lift
  has `v2 = k` exactly. Hence `D_{r*} ≠ D_{r*+2^k}` and the `r*` column is nonzero.
- Exactly one nonzero column ⟹ `rank(S_odd) = 1`. ∎

## Status of verification (numerical, this session)

- Structural column-builder = dense `S_odd`, max|diff| = 0.0 (k=9,10,11); every non-r* column
  exactly zero; full pipeline rank=1 with correct r* for k=3..18.
- Independent from-scratch exact integer rank = 1 (sympy) for k=3..11.
- Cyclic-shift identity 0 violations k=4..12; r* the unique exceptional residue (brute force k<=20,
  modular inverse k<=24).

## What this does NOT do

This is a structural fact about one block of the (unproven-uniform-gap) transfer operator. It does
not establish the spectral gap (Section 2 of SESSION_39 synthesis - that inductive argument is
refuted) and does not touch the deterministic-orbit transfer principle. It is a clean, real,
formalisable lemma - the one genuinely new theorem of the session.
```

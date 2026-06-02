# Lemma B: the collision bound, unconditional for all k (2026-06-01)

> **RETRACTION note (2026-06-02).** Lemma B is a correct ingredient of the operator spectral-gap proof.
> The downstream conclusion "certificate => no non-trivial cycles" is **false and withdrawn** (the
> `3x-1` control has the same certificate and known cycles). See [CYCLE_CLAIM_REFUTED.md](CYCLE_CLAIM_REFUTED.md).
> The collision bound itself stands; read any "eliminating cycles" phrasing below as withdrawn.

Lemma B is the "diagonal + lower back-flow" half of the uniform spectral-gap certificate
for the Syracuse transfer operator (see
[STEP4_BLOCK_FORMULA_FOUNDATION.md](STEP4_BLOCK_FORMULA_FOUNDATION.md)).

**Honest scope.** This document closes the *combinatorial core* of Lemma B - the bound
`coll(k) <= 3*2^k` - **unconditionally for all k**, replacing the previous finite-k verification.
The reduction from Lemma B to that bound uses three foundation facts R1, R2, R3: R3
(`||u||^2 = 1 - 2^{1-k}`) is unconditional; R1 (strict-upper `U_clean`) and R2 (defect `D` confined to
row `r*`, `rank(D) <= 1`) are also now **proved unconditionally for all k** - they reduce by elementary
algebra to the Coset-Uniformity lemma CU, proved by finite-group theory in
[HALFSHIFT_S4_LEMMA_A_PROOF.md](HALFSHIFT_S4_LEMMA_A_PROOF.md) Sections 1-2. (Earlier drafts marked
Half-Shift Invariance as "finite-k verified"; that is superseded - CU is the analytic proof. The only
fact that remains finite-k-verified is the defect *nonvanishing* / "rank exactly 1", which is not used:
every bound here is an upper bound that `D = 0` would satisfy trivially.) So Lemma B and its foundation
are unconditional; the collision bound below is its combinatorial heart.

The one load-bearing combinatorial gap (per-shell injectivity, "FACT 1") had a loose justification
in [attack3_provable_bound.py](attack3_provable_bound.py)
and was stated there for even k only (the `a = 1` model); the rigorous parity-correct argument and
its full ground-truth verification are in
[lemmaB_fact1_rigorous.py](lemmaB_fact1_rigorous.py).

## Statement

Let `Q_k` be the `(k-1) x (k-1)` matrix of 2-adic level-block operator norms
`Q_k[a,b] = ||P_a U_k P_b||_2`, `a, b in {0,...,k-2}`. Then

```
    || tril(Q_k) ||_2  <=  sqrt(3) * 2^{-k/2}     for all k,
```

where `tril` is the lower-triangular part (`a >= b`). The constant `sqrt(3)` is provable; the
sharp constant is `sqrt(31/12) = 1.6073` (even k, verified), well inside `sqrt(3) = 1.7321`.

## Reduction to a single collision count

Three facts from the block-formula foundation (all PROVED for all k there):

- **R1.** `U_k = U_clean + D`, with `U_clean` strictly upper-triangular in level
  (`P_a U_clean P_b = 0` for `a >= b`). Hence `tril(Q_k) = tril(Q_D)`: the entire lower
  triangle comes from the defect `D`.
- **R2.** `D = e_{r*} c^*` is rank-1, supported on the single row `r* = -3^{-1} mod 2^k`.
  So `P_a D P_b = (P_a e_{r*})(P_b c)^*` is rank-1 and
  `Q_D[a,b] = ||P_a D P_b||_2 = u_a v_b`, where `u_a = ||P_a e_{r*}||`, `v_b = ||P_b c||`.
- **R3.** `e_{r*}` is a state-basis unit vector, so `u_a^2 = d_a / N` with `d_a = 2^{k-2-a}`,
  `N = 2^{k-1}`. Summed over the levels `a = 0..k-2` (the Perron direction `xi=0` excluded):

```
    ||u||^2 = sum_{a=0}^{k-2} d_a / N = (N-1)/N = 1 - 2^{1-k} < 1.       [verified exactly]
```

  And `||v||^2 = sum_b v_b^2 = ||c||^2 - 1/N <= ||c||^2`: the level sum runs over `b = 0..k-2` and
  omits the Perron direction `xi = 0`, where `c` carries `|<chi_0, c>|^2 = 1/N` (since
  `sum_s c[s] = 1`). The omitted term only strengthens the bound.

Now bound the lower-triangular spectral norm by its Frobenius norm (a sum of fewer nonnegative
terms than the full outer product):

```
    || tril(Q_k) ||_2  =  || tril(Q_D) ||_2
                       <= || tril(Q_D) ||_F
                       <= || u v^* ||_F  =  ||u|| * ||v||
                       <= ||u|| * ||c||  <  ||c||.                        (since ||v|| <= ||c||, ||u|| < 1)
```

So Lemma B follows from a bound on `||c||` alone, and `||c||` is a collision count.

The defect covector is `c[s] = (1/2^k) #{m in Z/2^k : Syr(r* + m 2^k) = 2s+1 (mod 2^k)}`. Writing
`cf[t] = #{m : Syr(r*+m 2^k) = t}` (so `c = cf / 2^k` reindexed onto odd residues),

```
    ||c||^2 = sum_s c[s]^2 = (1/4^k) sum_t cf[t]^2 = coll / 4^k,
    coll := sum_t cf[t]^2.
```

Hence `||c|| = sqrt(coll) / 2^k`, and Lemma B reduces to:

```
    coll(k) <= 3 * 2^k    for all k.                                     [the whole content]
```

Then `||tril(Q_k)||_2 <= ||u|| ||c|| <= sqrt(coll)/2^k <= sqrt(3) 2^{-k/2}`.

## The collision bound  coll(k) <= 3 * 2^k

The fibre value is `Syr(r* + m 2^k) = oddpart(3(r*+m 2^k)+1) (mod 2^k)`. Since `r*` is the unique
odd residue with `v2(3r*+1) = 2 ceil(k/2)`, we have `3r*+1 = 2^{2 ceil(k/2)}`, so

```
    3(r* + m 2^k) + 1 = 2^{2 ceil(k/2)} + 3m 2^k = 2^k (a + 3m),
    a := 2^{2 ceil(k/2) - k} in {1, 2} :  a = 1 (k EVEN, 3r*+1 = 2^k),  a = 2 (k ODD, 3r*+1 = 2^{k+1}).
```

Hence `Syr(r* + m 2^k) = oddpart(a + 3m) (mod 2^k)`, where `x = a + 3m` ranges over the
**real-integer** arithmetic progression `{a, a+3, ..., a + 3(2^k - 1)}` of length `2^k`. This model
is cross-checked against the true Syracuse fibre `syr` for every k and both parities in
[lemmaB_fact1_rigorous.py](lemmaB_fact1_rigorous.py)
(column `=syr?` all `True`). (The even-k case, with `a = 1`, was the one established earlier in
[attack3_collision_proof.py](attack3_collision_proof.py);
the `a = 2` odd-k case is the corrected generalisation.)

Decompose the AP by 2-adic valuation into shells `A_j = { x = a+3m : v2(x) = j }`.

**Shell sizes.** Since `3` is a unit mod `2^{j+1}`, `x = a+3m mod 2^{j+1}` is equidistributed as
`m` runs over `2^k` consecutive AP terms, so `|A_j| = 2^{k-1-j}` for `j = 0..k-1`, plus a single
top atom `x ≡ 0 (mod 2^k)` (the unique `m` with `v2(x) >= k`). Total `= (2^k - 1) + 1 = 2^k`.

For `x in A_j`, `oddpart(x) = u = x / 2^j` is odd; the within-shell image is `u mod 2^k`. Let
`R_j = { (x/2^j) mod 2^k : x in A_j }`.

### FACT 1 (per-shell injectivity) -- rigorous

> The map `x -> (x / 2^j) mod 2^k` is injective on `A_j`, for every `j = 0..k-1`.
> Hence `|R_j| = |A_j| = 2^{k-1-j}` and `cf` restricted to shell `j` is 0/1.

**Proof.** Suppose `x = 2^j u`, `x' = 2^j u'` lie in `A_j` with `u ≡ u' (mod 2^k)` and `x ≠ x'`.

1. `u - u' = 2^k s` for some integer `s ≠ 0`, so `x - x' = 2^j(u - u') = 2^{k+j} s`.        (I)
2. `x ≡ x' ≡ a (mod 3)` (both are `a + 3m` terms, a single class mod 3), so
   `3 | (x - x') = 2^{k+j} s`; as `gcd(3, 2) = 1`, this forces `3 | s`.                       (II)
3. `x, x' in [a, a + 3*2^k - 3] ⊂ [lemmaB_fact1_rigorous.py](lemmaB_fact1_rigorous.py)
(columns `(I) (II) (R) inj` all `True`).

### Diagonal and cross terms

By FACT 1, `cf[t] = #{ j : t in R_j }` (the number of shells whose residue set contains `t`).
Therefore

```
    coll = sum_t cf[t]^2 = sum_t (#shells hitting t)^2 = sum_{j,j'} |R_j ∩ R_{j'}|
         = diag + 2 * cross,
    diag  = sum_j |R_j|,
    cross = sum_{j < j'} |R_j ∩ R_{j'}|.
```

- **Diagonal (exact).** `diag = sum_j |R_j| = (2^k - 1) + 1 = 2^k`. The `+1` is the top atom
  (`oddpart 0`, hit once). PROVED, verified `diag = 2^k` for all tested k.
- **Cross (bound).** For `j < j'`, `|R_j ∩ R_{j'}| <= |R_{j'}| = 2^{k-1-j'}`. Each `j'` has `j'`
  lower partners `j = 0..j'-1`, so

```
    cross <= sum_{j'=0}^{k-1} j' * 2^{k-1-j'}
          <= 2^{k-1} * sum_{j'>=0} j' 2^{-j'}  =  2^{k-1} * 2  =  2^k.
```

Hence `coll = diag + 2*cross <= 2^k + 2*2^k = 3 * 2^k`. **QED, all k.**

### Sharp constant (context, not needed)

`coll(k)/2^k -> 31/12 = 2.5833` for both parities, so `||c|| * 2^{k/2} -> sqrt(31/12) = 1.6073`.
The even-k exact count is `coll = (31*2^k + 8)/12` (confirmed `166, 662, 2646, 10582` at
k=6,8,10,12); the odd-k count approaches the same ratio from below (`80, 328, 1320, 5288, 21160`
at k=5,7,9,11,13, all matching the true Syracuse fibre). The provable `sqrt(3) = 1.7321` bound has
comfortable margin over `sqrt(31/12)` at every k.

## What this closes

The collision bound `coll <= 3*2^k` is now UNCONDITIONAL for all k (both parities), from the
parity-correct shell decomposition: FACT 1 (per-shell injectivity, the mod-3 + range argument),
`diag = 2^k` exact, `cross <= 2^k`. Combined with R3 (`||u||^2 = 1 - 2^{1-k}`) and the Frobenius
bound, this gives `||u|| ||c|| <= sqrt(3) 2^{-k/2}`.

The remaining step to `||tril(Q_k)||_2 <= sqrt(3) 2^{-k/2}` is the identity `tril(Q_k) = tril(Q_D)`
and `Q_D[a,b] = u_a v_b`, i.e. R1 (strict-upper `U_clean`) and R2 (defect `D` confined to the row
`r*`, hence `rank(D) <= 1`). **Update (2026-06-02): R1 and R2 are now proved unconditionally for all
k** - they reduce by elementary algebra to the Coset-Uniformity lemma CU, which is proved (finite-group
theory, no finite-k input) in [HALFSHIFT_S4_LEMMA_A_PROOF.md](HALFSHIFT_S4_LEMMA_A_PROOF.md) Sections
1-2. (The Frobenius bound uses only `rank(D) <= 1` / the `D = e_{r*} c^*` *form*, an upper bound; the
separate "rank exactly 1" / defect-nonvanishing fact is verified but not analytically settled for even
k, and is NOT used - a zero defect would only make `||tril Q_k|| = 0`.) So `||tril(Q_k)||_2 -> 0` at
rate `2^{-k/2}` holds for all k unconditionally.

**Lemma A** (the upper cascade `||P_a U_k P_b|| = 2^{-(b-a)/2}` via the within-level isometry
`B*B = 2^{-(b-a)} I`) is likewise now proved for all k (CU + Shell-Bijection + S4, same document,
Section 4). With Lemmas A, B and the per-level Lemma C ([LEMMA_C_PROOF.md](LEMMA_C_PROOF.md)), the
certificate row-sum is `< G_up + 2^{-3/2} = 0.9005 < 1` for all k ([UFULL_ASSEMBLY_PROOF.md](UFULL_ASSEMBLY_PROOF.md)),
eliminating non-trivial Syracuse cycles. The remaining program-level obligations are the explicit
write-up of the standard spectral reduction (`cert < 1 => |lambda_2(U_k)| < 1`) and the Lean
formalisation - not any finite-k crux.

## Verification

```
python lemmaB_fact1_rigorous.py     # parity-correct FACT 1 (I/II/R) + diag/cross/coll + Frobenius
                                    #   chain, with =syr? ground-truth cross-check, both parities
python adv_tril_sep_correct.py      # ||tril(Q_D)|| matches dense build_T (reldiff ~1e-11), chain k=6..24
python attack3_provable_bound.py    # original shell argument (NOTE: even-k a=1 model only; odd k
                                    #   understated - superseded by lemmaB_fact1_rigorous.py)
```

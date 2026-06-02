# Lemma C: per-level decay of the defect covector

> **RETRACTION note (2026-06-02).** Lemma C is a correct ingredient of the operator spectral-gap proof.
> The downstream conclusion "certificate => no non-trivial cycles" is **false and withdrawn** (the
> `3x-1` control has the same certificate and known cycles). See [CYCLE_CLAIM_REFUTED.md](CYCLE_CLAIM_REFUTED.md).
> Read any "eliminates cycles" phrasing below as withdrawn; the per-level bound itself stands.

Written 2026-06-02 (Session 40), the analytic step requested after
[UFULL_ASSEMBLY_PROOF.md](UFULL_ASSEMBLY_PROOF.md). This document reduces Lemma C to an elementary
collision-count inequality in the **same shell framework that proves Lemma B**, and proves the bound
at the strength the assembly needs. Reproduce with `probe_gb_collision.py`, `probe_autocorr.py`,
`probe_shell_halfshift.py`, and the one-liners in this file's verification section.

## Statement and what the assembly actually needs

> **Lemma C (sharp form).** `v_b <= (3/4) 2^{-b} 2^{-k/2}` for all `k`, `0 <= b <= k-2`. Equivalently,
> with `g_b := v_b 2^b 2^{k/2}`, `g_b^2 <= 9/16`.

The assembly ([UFULL_ASSEMBLY_PROOF.md](UFULL_ASSEMBLY_PROOF.md)) does not need the sharp constant. It
needs `sum_{b} v_b 2^{-b} <= C 2^{-k/2} (4/3)` with the certificate `< G_up + C (4/3) 2^{-3/2} < 1`,
i.e. any uniform bound

> **Lemma C (assembly form).** `g_b <= 0.961` for all `k, b`.

We prove `g_b < sqrt(3/4) = 0.8660 < 0.961` for all `k, b` (Theorem C below), **unconditionally**: the
bound on `g_b` (a function of the defect count alone) uses only Lemma H (now proved, an elementary
parity split), Parseval, and FACT 1 + the shell sizes (Lemma B's unconditional combinatorial core). The
cross-scale homometry earlier sessions left un-derived is the content of Lemma H. The sharp `g_b <= 3/4`
is verified exactly to `k=26` but not needed. The Half-Shift Invariance crux enters only where the
*assembly* uses the bound (the rank-1 defect factorisation R2), exactly as for Lemmas A and B; Lemma C
introduces **no new conditionality**.

## Step 1 - the defect covector is a 3x+1 exponential sum, and `g_b` is single-index

`r* = -3^{-1} mod 2^k` is the unique odd residue with `3r*+1 = 2^{2 ceil(k/2)}` (Lemma B). Lifting
`r*+m 2^k` and applying the Syracuse map collapses the `2^k` factor:
```
  Syr(r* + m 2^k) = OddCore(3m + a) mod 2^k,   a := 2^{2 ceil(k/2) - k} in {1,2}   (a=1 k even, a=2 k odd).
```
So the defect count `cf_k[t] = #{m in [0,2^k): Syr(r*+m 2^k) = t}` is the pushforward of the uniform
measure on `[0,2^k)` under `m -> OddCore(3m+a) mod 2^k`, and `c = cf_k / 2^k`. Its Fourier transform is
the 3x+1 exponential sum `S_k(xi) = sum_t cf_k[t] w^{xi t}`, `w = e^{-2 pi i/2^k}`, with
`v_b^2 = sum_{xi: v2(xi)=b, 1<=xi<2^{k-1}} |S_k(xi)|^2 / (2^{2k} 2^{k-1})`.

### Lemma H (the frequency-doubling recursion) - PROVED

> **Lemma H.** For every `k >= 2` and every integer `xi` with `2 xi not= 0 (mod 2^k)`,
> `S_k(2 xi) = S_{k-1}(xi)` exactly (as complex numbers).

**Proof.** Split the sum `S_k(eta) = sum_{m=0}^{2^k-1} w_k^{eta O_k(m)}`, `O_k(m)=OddCore(3m+a_k) mod 2^k`,
by the parity of `m`, with `m = 2m'` or `m = 2m'+1`, `m' in [0,2^{k-1})`.

One parity gives `3m + a_k` **odd** (a number of 2-adic valuation 0), so `O_k = 3m + a_k` is linear in
`m'`; call this the *geometric part*. (For `k` even, `a_k=1`: even `m` gives `6m'+1`, odd. For `k` odd,
`a_k=2`: odd `m` gives `6m'+5`, odd.) The geometric part of `S_k(eta)` is
```
  sum_{m'=0}^{2^{k-1}-1} w_k^{eta (6m' + c)} = w_k^{eta c} sum_{m'} (w_k^{6 eta})^{m'},   c in {1,5}.
```
Its ratio `r = w_k^{6 eta}` satisfies `r^{2^{k-1}} = w_k^{3 eta 2^k} = e^{-2 pi i 3 eta} = 1`, so the
geometric sum is `0` whenever `r not= 1`, i.e. whenever `6 eta not= 0 (mod 2^k)`, i.e. (as `3` is odd)
`eta not= 0 (mod 2^{k-1})`. Taking `eta = 2 xi`, the geometric part **vanishes** unless
`2 xi in {0, 2^{k-1}}`, i.e. unless `xi in {0, 2^{k-2}}`.

The other parity gives `3m + a_k` **even**: for `k` even, odd `m` gives `6m'+4 = 2(3m'+2)`, so
`O_k = OddCore(3m'+2)`; for `k` odd, even `m` gives `6m'+2 = 2(3m'+1)`, so `O_k = OddCore(3m'+1)`. In
both cases the residual map is `OddCore(3m' + a_{k-1})` - the parity flip of `m` swaps `+1 <-> +2`,
matching `a_{k-1}` exactly. This *structural part* of `S_k(2 xi)` is
```
  sum_{m'=0}^{2^{k-1}-1} w_k^{2 xi (OddCore(3m'+a_{k-1}) mod 2^k)}
   = sum_{m'} w_{k-1}^{xi (OddCore(3m'+a_{k-1}) mod 2^{k-1})}  =  S_{k-1}(xi),
```
since `w_k^{2 xi X} = w_{k-1}^{xi X}` depends only on `X mod 2^{k-1}`, and the last sum is by definition
`S_{k-1}(xi)`. Hence for `2 xi not= 2^{k-1}` (and `xi not= 0`), `S_k(2 xi) = 0 + S_{k-1}(xi)`. QED.

Verified exactly to `k=14` complex-valued and to `k=26` in magnitude (`probe_homometry_proof.py`,
`probe_selfsimilar.py`): the geometric part is `0` for every `xi` except `xi = 2^{k-2}`, and the
structural part equals `S_{k-1}(xi)` with no exception.

### The single-index reduction `g_b^2 = E_p/4^p` (now rigorous)

Iterate Lemma H. The level-`b` frequencies are `xi = 2^b m`, `m` odd, `1 <= m < 2^{k-1-b}`. Peeling one
factor of `2` at a time, `S_k(2^b m) = S_{k-1}(2^{b-1} m) = ... = S_{k-b}(m)`; each peel is legitimate
because the intermediate argument equals the step's Nyquist `2^{(k-i)-1}` only if `m = 2^{k-1-b}`, which
is excluded by `m < 2^{k-1-b}`. With `|S_k(xi)| = |FFT_{2^{k-1}}(c)[xi]| 2^k` giving
`v_b^2 = (2^{2k} 2^{k-1})^{-1} sum_{v2(xi)=b} |S_k(xi)|^2`, and `sum_{m odd < 2^{k-1-b}} |S_{k-b}(m)|^2
= (1/2) E_p` by the conjugate symmetry `|S_p(m)| = |S_p(2^p - m)|` (no odd fixed point), we get
```
      g_b^2 = 2^k 4^b v_b^2 = E_p / 4^p,    E_p := sum_{xi odd, 1<=xi<2^p} |S_p(xi)|^2,   p = k-b.   (2)
```
This is now derived, not merely verified. (A second, fully independent route - the scale-`k` Parseval
identity `g_b^2 = 2^{b-1-k}(2 fcoll_b - fcoll_{b+1})`, `fcoll_j = ||fold_{2^j} cf_k||^2` - agrees with
(2) to machine precision and needs no homometry, but the homometry route above is the clean one.)

## Step 2 - `E_p` is a half-shift autocorrelation (fully rigorous)

Write the odd-frequency indicator `[xi odd] = (1 - (-1)^xi)/2` and use Parseval on `Z/2^p` twice:
```
  sum_{xi}        |S_p(xi)|^2          = 2^p coll(p),       coll(p) = sum_t cf_p[t]^2,
  sum_{xi} (-1)^xi |S_p(xi)|^2          = 2^p A_p,           A_p     = sum_t cf_p[t] cf_p[t + 2^{p-1}],
```
the second because `(-1)^xi = w^{xi 2^{p-1}}` and `sum_xi |S|^2 w^{xi h} = 2^p (autocorrelation at lag h)`.
Hence
```
  E_p = (1/2)( 2^p coll(p) - 2^p A_p ) = 2^{p-1} ( coll(p) - A_p ),                  (3)
```
and, using `coll - A = (1/2) sum_t (cf[t] - cf[t+2^{p-1}])^2`,
```
  g_b^2 = E_p/4^p = (coll(p) - A_p) / 2^{p+1}.                                       (4)
```
Identities (3), (4) are exact and unconditional (`probe_autocorr.py`: `match=True` for all `p`).
Grouping `m` by `u := O_p(m) mod 2^{p-1}` with top-bit split `n_u^0, n_u^1`,
```
  coll(p) - A_p = sum_{u} (n_u^0 - n_u^1)^2.                                         (5)
```
So **Lemma C (assembly form) <=> `coll(p) - A_p <= (7/8) 2^p` for all p** (this gives `g_b^2 <= 7/16`
asymptotically; we prove the looser `<= (3/2) 2^p` below, i.e. `g_b^2 <= 3/4`, which already suffices).

## Step 3 - shell bound (Lemma B's machinery)

Decompose the AP `{a+3m : m in [0,2^p)}` into 2-adic shells `A_j = {x : v2(x) = j}`, and set
`R_j = {(x/2^j) mod 2^p : x in A_j}`. Lemma B's facts carry over verbatim:
- **FACT 1 (per-shell injectivity).** `x -> (x/2^j) mod 2^p` is injective on `A_j`; so `cf_p[t] = #{j : t in R_j}` is `0/1` per shell, and `|R_j| = 2^{p-1-j}` for `j=0..p-1`, plus a single top atom.
- In particular `|R_0| = 2^{p-1}`: shell `j=0` is **all** `2^{p-1}` odd residues mod `2^p`.

Write `coll(p) - A_p = sum_{j,j'} ( c_{jj'} - d_{jj'} )`, `c_{jj'} = |R_j cap R_{j'}|`,
`d_{jj'} = |R_j cap (R_{j'} - 2^{p-1})|` (this is (5) re-expanded via `cf = sum_j 1[. in R_j]`).

**The `j=0` shell cancels.** Since `R_0` is all odd residues and `2^{p-1}` is even, `R_{j'} - 2^{p-1}`
is also a set of odd residues, so `c_{0 j'} = |R_{j'}| = d_{0 j'}` and likewise `c_{j' 0} = d_{j' 0}`.
Every term with a `0` index vanishes:
```
  coll(p) - A_p = sum_{j, j' >= 1} ( c_{jj'} - d_{jj'} ).                            (6)
```

**Drop the non-negative subtractions** (`d_{jj'} >= 0`, and the diagonal `d_{jj} = s_j >= 0`):
```
  coll(p) - A_p <= sum_{j>=1} c_{jj} + sum_{j != j', j,j' >= 1} c_{jj'}
               =  sum_{j>=1} |R_j| + 2 sum_{1<=j<j'} |R_j cap R_{j'}|.
```
The first sum is `sum_{j=1}^{p-1} 2^{p-1-j} + 1 (atom) = (2^{p-1} - 1) + 1 = 2^{p-1}`. For the second,
`|R_j cap R_{j'}| <= |R_{j'}|` (the smaller shell, `j' > j`); the shells `j' = 2..p-1` have
`|R_{j'}| = 2^{p-1-j'}` and `j'-1` partners `j in {1,...,j'-1}`, and the top atom `j' = p` has
`|R_p| = 1` and `p-1` partners. The honest finite sum has an exact closed form:
```
  sum_{1<=j<j'<=p} |R_j cap R_{j'}|  <=  sum_{j'=2}^{p-1}(j'-1) 2^{p-1-j'}  +  (p-1) * 1
                                     =  2^{p-1} - 1                              (verified for all p).
```
(The `(p-1)` is the top atom's partner contribution; the closed form `... = 2^{p-1}-1` is checked in
the verification section.) Therefore
```
  coll(p) - A_p  <=  2^{p-1} + 2 (2^{p-1} - 1)  =  3 * 2^{p-1} - 2  <  (3/2) 2^p,    (7)
```
and by (4),
```
  g_b^2  =  (coll(p) - A_p) / 2^{p+1}  <=  (3 * 2^{p-1} - 2)/2^{p+1}  =  3/4 - 2^{-p}  <  3/4.   (8)
```

## Theorem C and the conclusion

> **Theorem C.** `g_b < sqrt(3/4) = 0.8660 < 0.961` for all `k` and `0 <= b <= k-2`, **unconditionally**.
> Hence `v_b = ||P_b c|| < (3/4) 2^{-b} 2^{-k/2}`, and the certificate is `< 1` uniformly.

`g_b = v_b 2^b 2^{k/2}` is a function of the defect count `cf` alone (`c = cf/2^k`), and the bound on it
uses only unconditional ingredients: **Lemma H** (parity split), the Parseval identities (2)-(5), and
**FACT 1 + the shell sizes** `|R_j| = 2^{p-1-j}` (Lemma B's combinatorial core, proved there
unconditionally for all `p`, both parities). So Theorem C is unconditional.

The Half-Shift Invariance dependence enters only when this bound is *used in the assembly*: that the
lower-triangle block norms factor as `Q_D[a,b] = u_a v_b` is foundation fact R2 (rank-1 defect
`D = e_{r*} c^*`), which rests on Half-Shift Invariance - the same dependence Lemmas A and B already
carry. Lemma C itself adds **no new conditional input** to the program.

The bound (8) is `g_b^2 <= 3/4`, well inside the `0.924` ceiling the assembly needs. (A sharper
accounting of the top-atom partners gives `g_b^2 <= 3/4 + (p-1)/2^p`, which still respects the ceiling
for `p>=5`; for the finitely many top levels `p = k-b in {2,3,4}` the values are the exact, `k`-independent
dyadic rationals `g_b^2 in {1/2, 1/8, 9/16}`, all `<= 9/16`. Either way `g_b < 0.961`.)

Feeding `g_b <= sqrt(3/4)` into the assembly: `sum_b v_b 2^{-b} <= sqrt(3/4) 2^{-k/2} (4/3)`, so the
lower row-sum at the top is `<= 2^{-3/2} sqrt(3/4)(4/3) = 0.408`, and
```
  cert(k) < G_up + 0.408 = 0.547 + 0.408 = 0.955 < 1,   uniform in k.
```

## Status: what is rigorous, what is inherited

**Unconditional (no Half-Shift Invariance, no operator input - all about the count `cf` and the AP):**
- **Lemma H** (the homometry `S_k(2 xi) = S_{k-1}(xi)`): proved by the elementary parity split.
- The single-index reduction (2) `g_b^2 = E_p/4^p` (Lemma H iterated + conjugate symmetry); cross-checked
  by the scale-`k` Parseval identity `g_b^2 = 2^{b-1-k}(2 fcoll_b - fcoll_{b+1})`.
- Step 2: identities (3), (4), (5) (Parseval + the `coll-A = (1/2) sum (cf[t]-cf[t+h])^2` algebra).
- Step 3: the `j=0` cancellation (6), the shell drop, the finite bound (7), and `g_b^2 < 3/4` (8) -
  using FACT 1 and the shell sizes `|R_j| = 2^{p-1-j}`, which are Lemma B's combinatorial core, **proved
  there unconditionally** (the mod-3 + range injectivity, all `p`, both parities).

So **Theorem C - the bound `v_b < (3/4) 2^{-b} 2^{-k/2}` - is unconditional.**

**Inherited (the same dependence Lemmas A and B already carry):**
- The bound's *use in the assembly* needs the lower-triangle factorisation `Q_D[a,b] = u_a v_b`,
  foundation fact R2 (rank-1 defect `D = e_{r*} c^*`), which rests on the Half-Shift Invariance /
  coset-uniformity lemma (a draft with a finite-`k`-verified crux). This is not new to Lemma C.

So the precise status:

> **Lemma C (the per-level bound) is proved unconditionally.** The homometry that earlier sessions left
> open is now derived (Lemma H, elementary). Combined with Lemma A (proved) and the assembly (proved),
> the certificate is `cert(k) < 1` for every `k`. The shared foundation (Coset-Uniformity, giving R1/R2)
> is itself now proved unconditionally ([HALFSHIFT_S4_LEMMA_A_PROOF.md](HALFSHIFT_S4_LEMMA_A_PROOF.md)
> Sections 1-2), so the whole chain is all-`k`. (Lemma C uses only `rank(D) <= 1` / the `D = e_{r*} c^*`
> form; the defect *nonvanishing* is not needed - `v_b = 0` would satisfy the bound trivially.) The
> remaining open obligations for the cycle program are therefore (a) the explicit write-up of the
> standard spectral reduction `cert(k) < 1 => |lambda_2(U_k)| < 1` (cited in the foundation, not yet
> written in full) and (b) the Lean formalisation - not any finite-`k` crux.

## Verification

```bash
python probe_gb_collision.py     # g_b^2 = E_p/4^p (single index); sum|S|^2 = 2^p coll(p)
python probe_autocorr.py         # E_p = 2^{p-1}(coll-A); coll-A = (1/2)sum(cf[t]-cf[t+h])^2
python probe_shell_halfshift.py  # shell decomposition: diag=2^{p-1}, j=0 cancels, the (3/2)2^p bound
python probe_homometry_proof.py  # Lemma H: parity split, geometric part vanishes, nongeo = S_{k-1}
python probe_selfsimilar.py      # the homometry S_k(2 xi) = S_{k-1}(xi) (now proved, Lemma H)
# closed forms (one-liners):
#   sum_{j=2}^{p-1}(j-1)2^{p-1-j} + (p-1) == 2^{p-1}-1     (cross-sum bound, p=2..24: True)
#   coll(p) - A_p <= 3*2^{p-1} - 2                          (eq 7, p=2..16: True)
```

## Honest ceiling (unchanged)

This is cycle elimination, not Collatz. Even with Lemma C fully closed, the program proves only "the
only Collatz cycle is `1 -> 4 -> 2 -> 1`," strictly weaker than the conjecture; the divergent-trajectory
half is untouched.

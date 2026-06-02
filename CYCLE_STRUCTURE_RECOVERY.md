# Can cycle-specific structure be recovered? (investigation)

Written 2026-06-02 (Session 40), following the retraction in [CYCLE_CLAIM_REFUTED.md](CYCLE_CLAIM_REFUTED.md).
Question: after the spectral gap was shown blind to cycles, is there *any* refinement of the mod-`2^k`
transfer-operator framework that recovers cycle-specific information? Tested empirically on
`3x+1` (no nontrivial cycle) vs `3x-1` (cycles `{5,7}`, `{17,25,37,41,55,61,91}`), where ground truth
is known. Reproduce: `python probe_cycle_recovery.py`.

## What does NOT recover cycles (verified cycle-blind)

| Candidate invariant | 3x+1 vs 3x-1 | Verdict |
|---|---|---|
| Spectral gap `\|lambda_2\|` | 0.25 vs 0.29 (both < 1) | blind (already known) |
| Eigenvalues near `L`-th roots of unity | nearest non-Perron eigenvalue ~0.72 from any root, **both** maps; the `{5,7}` 2-cycle yields **no** eigenvalue near `-1` | blind |
| Traces `tr(T_k^n)`, `n=1..6` | all `~1.000` for both (Perron-dominated; cycle weight is `O(2^{-k})`) | blind |
| Row-sum certificate `cert(k)` | 0.634 vs 0.606 (both < 1) | blind |

The reason is structural and was the content of the retraction: `T_k` averages each residue over all
`2^k` lifts (coset-uniformity), giving each true cycle edge weight `2^{-k}` and scattering the rest.
Every *spectral* functional of `T_k` is therefore insensitive to cycles, confirmed across the gap,
the full spectrum, roots of unity, and traces.

## Where cycle structure actually lives (and how it touches the framework)

A cycle `a_1 -> ... -> a_L -> a_1` with valuations `v_i = v2(3 a_i + s)` (`s = +-1`) satisfies the
**cycle equation**
```
    a_1 (2^V - 3^L) = s * sum_{i=1}^{L} 3^{L-i} 2^{v_1 + ... + v_{i-1}},   V = v_1 + ... + v_L.
```
Measured for the real cycles:

| cycle | L | V | V/L | `2^V - 3^L` |
|---|---|---|---|---|
| 3x-1 `{5,7}` | 2 | 3 | 1.50 | `-1` |
| 3x-1 7-cycle | 7 | 11 | 1.571 | `-139` |
| 3x+1 `{1}` (trivial) | 1 | 2 | 2.00 | `+1` |
| 3x-1 `{1}` | 1 | 1 | 1.00 | `-1` |

Two facts this makes precise:

1. **The discriminating quantity is `2^V - 3^L`, and its SIGN is the `+1` vs `-1` difference.** For
   `s=+1` a cycle needs `2^V - 3^L > 0`; for `s=-1`, `< 0`. Operator norms / singular values / `|Q[a,b]|`
   are sign-blind (they depend on `|w^{...}|`-type magnitudes), so they cannot see this. This is exactly
   why the identical-looking certificate certifies both maps.

2. **The valuations `v_i` ARE the framework's 2-adic levels.** `v_i = v2(3 a_i + s)` is precisely the
   level `v(r)` that indexes the operator's shells and the coset-uniformity lemma. So the framework
   shares *ingredients* with the cycle equation - but the cycle question is the Diophantine balance
   `2^V approx 3^L` (a linear form in logs `V log 2 - L log 3`), not a spectral quantity. A cycle is a
   rare sustained run of below-average valuations (`V/L approx log2(3) = 1.585`, vs the mean `E[v] = 2`);
   the gap measures the *average* drift (the same `E[v]=2` for both maps), not the rare balanced runs.

## The intrinsic tension

- **Cycle-sensitive operator = the deterministic map.** The only construction here that separates the
  two maps (probe `(C)`: 1 cycle vs 3) is the deterministic partial map `D_k(n) = Syr(n)` on
  `n < 2^k` - i.e. literally iterating Syracuse. Its periodic orbits are the cycles; "computing" them
  is brute-force search. The Ruelle/zeta-function transfer operator (sum over inverse branches,
  weighted by contraction) is the principled version - its traces count periodic orbits - but it lives
  on an infinite space and proving "no nontrivial cycle" via its zeta function is as hard as the
  problem itself.
- **Tractable operator = cycle-blind.** The averaged `T_k` is analytically tractable *because* the
  coset-averaging gives coset-uniformity and the clean spectral gap. That same averaging erases cycles.

So tractability and cycle-sensitivity are, in this framework, mutually exclusive: the property that
makes the spectral gap provable is the property that makes it blind to cycles.

## Conclusion

**No cycle-specific structure is recoverable from the spectral side of the transfer-operator framework.**
What the framework *does* supply cleanly - the 2-adic valuation / shell structure (coset-uniformity,
the `r*` defect, the level decomposition) - are the `v_i` ingredients of the cycle equation, so the
framework is not irrelevant to cycles; but the actual obstruction to a cycle is the Diophantine
inequality on `2^V - 3^L` (linear forms in logs, Baker-type bounds), which is a different toolset. For
`3x+1` that route is the source of the known results (no cycles below ~`10^{20}`; any cycle must have
astronomically large length), and it is sign-sensitive in exactly the way the spectral gap is not.

If a transfer-operator route to cycle elimination exists at all, it would have to be built on a
cycle-sensitive (deterministic / Ruelle-trace) operator and would need a genuinely new idea to bound
periodic-orbit traces - not the averaged operator whose gap is proved here.

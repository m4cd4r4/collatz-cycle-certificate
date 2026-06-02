# A uniform spectral gap for the Syracuse transfer operator

> **RETRACTION (2026-06-02, read first).** This repository previously claimed that the spectral-gap
> certificate **eliminates non-trivial Collatz cycles**. **That claim is withdrawn - it is false.** A
> uniform spectral gap of the transfer operator does NOT imply the absence of cycles. Decisive control:
> the `3x-1` map has known non-trivial cycles (`{5,7}`, `{17,25,...}`) yet its transfer operator passes
> the identical certificate (`cert ~ 0.606 < 1`, `|lambda_2| ~ 0.29`) even more strongly than `3x+1`.
> The averaging over `2^k` lifts that defines the operator washes out the deterministic orbit
> structure, so cycles are invisible to the spectrum (the same way the doubling map is mixing yet has
> dense periodic points). Full account: [CYCLE_CLAIM_REFUTED.md](CYCLE_CLAIM_REFUTED.md).

What remains correct and proved is a **uniform spectral gap for the Syracuse (3n+1) mod-`2^k` transfer
operator** `T_k`: the three lemmas (A, B, C), their shared Coset-Uniformity foundation, and the row-sum
assembly are proved for every scale `k` by elementary 2-adic / finite-group arguments, giving
`cert(k) < 0.9005 < 1` and hence `|lambda_2(T_k)| < 1` uniformly. This is a real result about the
operator. It simply does **not** carry the cycle-elimination corollary that was claimed, and it is not
a proof of (any part of) the Collatz conjecture.

> **Honest scope.** Cycle elimination is what was *attempted* and is now retracted. The divergent-
> trajectory half of Collatz was never addressed. The surviving content is the operator spectral gap.

![certificate vs k](figures/fig4_certificate.png)

*The certificate value stays flat below 1 as the scale `k` grows, giving a uniform spectral gap. (Note:
the `3x-1` operator, which has cycles, produces the same below-1 curve - so this does not forbid cycles;
see the retraction. Computed from the real operator; reproduce with `python generate_figures.py`.)*

---

## In plain terms

**The Collatz game.** Pick a whole number. If it is even, halve it. If it is odd, triple it and add
one. Repeat. The conjecture (open since the 1930s) says you always eventually reach 1. Two ways it
*could* fail: a number could **loop forever** in a cycle that never hits 1, or it could **grow to
infinity**. This project set out to attack the first failure mode (ruling out loops) - and the central
idea, it turns out, **does not work**. What survives is a clean side-result about a matrix.

**The trick that was tried.** Instead of following one number, follow the *cloud* of where numbers
land, and write that as a big table of numbers (a "transfer operator"). The hope was: if the operator's
second-largest characteristic size (eigenvalue) is **below 1** - a "spectral gap" - then the cloud
mixes and there is no room for a hidden loop. We proved that gap, with room to spare, at every scale.

**Why the idea fails.** A spectral gap measures mixing of the *averaged* cloud, and averaging throws
away the exact step-by-step orbit - which is the only thing a loop lives in. The clean test: the very
similar `3x-1` map (triple and subtract one) **does** have loops, like `5 -> 7 -> 5`. Yet its operator
has the **same** gap and passes the **same** certificate. So a gap below 1 cannot be what rules out
loops - if it were, it would wrongly "rule out" the `3x-1` loops that plainly exist. (Same lesson as
the doubling map `x -> 2x`, which mixes perfectly yet is full of loops.) The cycle claim is **withdrawn**;
details in [CYCLE_CLAIM_REFUTED.md](CYCLE_CLAIM_REFUTED.md).

![spectrum](figures/fig2_spectrum.png)

*The operator's eigenvalues. The top one is 1; the rest sit well below it. This gap is real and proved -
but it does not rule out cycles (the cyclic `3x-1` operator looks identical).*

**What we actually proved.** That spectral gap, for every scale `k`. The matrix splits into blocks by
how divisible by two a number is; the crux was showing each block is a clean, rigid rescaling
(an *isometry*). People expected deep "Gauss sum" machinery; it comes down to the single fact that
**3 is an odd number** (a unit modulo any power of two), plus bookkeeping. Three lemmas (A, B, C), their
shared foundation (coset-uniformity), and the final "below 1" assembly are all proved for every scale.

![block structure](figures/fig3_incidence.png)

*One block, drawn as its nonzero pattern: exactly one mark per row. That single combinatorial fact is
the entire reason the block is a rigid rescaling.*

**Where it stands.** A correct, elementary, all-scales proof of a uniform spectral gap for the Syracuse
transfer operator (`< 0.9005`, room to spare). The cycle-elimination conclusion that motivated it is
withdrawn. This is **not** a proof of any part of Collatz.

---

## The result, precisely

Let `U_k` be the Syracuse transfer operator on odd residues mod `2^k`, in the character (Fourier)
basis, decomposed by 2-adic valuation level. Write `Q_k[a,b] = ||P_a U_k P_b||_2` for the operator
norm of the level-`(a,b)` block. The **certificate** is

```
    rho(U_k) restricted to the non-Perron part  <=  rho(Q_k)  <=  max_a sum_b Q_k[a,b] * 2^{a-b}  <  1 ,
```

uniform in `k`, giving a uniform spectral gap `|lambda_2(U_k)| < 1`. **This gap does not eliminate
cycles** - it holds equally for the `3x-1` map, which has cycles ([CYCLE_CLAIM_REFUTED.md](CYCLE_CLAIM_REFUTED.md)).
The certificate is a true statement about the operator; the cycle inference drawn from it was wrong.

![Q heatmap](figures/fig1_Q_heatmap.png)

*The block-norm matrix `Q` (log scale). The bright upper triangle is the cascade `Q[a,b] = 2^{-(b-a)/2}`;
the dark lower triangle is the rank-1 `r*` defect, which vanishes like `2^{-k/2}`.*

The certificate splits into three lemmas (A, B, C), all **proved for all `k`** (elementary; Lean
formalisation pending), plus the row-sum assembly, also proved.

### Lemma A - upper cascade (within-level isometry)

For `0 <= a < b <= k-2`, `||P_a U_clean P_b||_2 = 2^{-(b-a)/2}` exactly, because the block `B` satisfies
`B*B = 2^{-(b-a)} I`. This reduces - by elementary algebra - to **two uniform-fibre lemmas sharing one
engine** ("3 is a unit mod `2^k` + a multiplication map has uniform fibres on a cyclic 2-group"):

- **CU** (coset-uniformity): also discharges Half-Shift Invariance and the rank-1 `S_odd` structure.
- **SB** (shell-bijection `r -> (3r+1)/2^j mod 2^{k-j}`): the instance the cascade uses.

plus **S4**, a one-line finite geometric series (`Sodd` vanishes iff `k-m <= v2(alpha) <= k-2`), which
was long mistaken for an analytic obstruction. Then `B*B = 2^{-d} I` follows from one nonzero per row
(disjoint column supports) + constant modulus. Verified exact to `k=16` (rational arithmetic, zero
error on and off diagonal), no parity dependence. Full proof: [HALFSHIFT_S4_LEMMA_A_PROOF.md](HALFSHIFT_S4_LEMMA_A_PROOF.md).

### Lemma B - lower back-flow (rank-1 defect)

`||tril(Q_k)||_2 <= sqrt(3) * 2^{-k/2} -> 0`. The lower triangle is entirely the rank-1 defect at the
exceptional residue `r* = -3^{-1} mod 2^k`; the bound reduces to a collision count `coll(k) <= 3*2^k`,
proved unconditionally for all `k` by a 2-adic shell decomposition (per-shell injectivity via a
`mod 3` + range argument). Sharp constant `coll/2^k -> 31/12`. Full proof: [LEMMA_B_PROOF.md](LEMMA_B_PROOF.md).

### Lemma C - per-level decay of the defect covector

`v_b := ||P_b c||_2 <= (3/4) * 2^{-b} * 2^{-k/2}` for all `0 <= b <= k-2`, all `k`. The constant `3/4`
is sharp (equality at `b = k-4`), verified exactly to `k=26`, with a `k`-independent boundary profile.
It reduces to a periodization-excess bound on the defect covector `c` (the partial Gauss sum at `r*`):
with `h_j = N ||fold_{2^j}(c)||^2`, the bound is `2^b (2 h_b - h_{b+1}) <= 9/16`. The `L^2` mass alone
(Lemma B) is **not** enough - it gives only `0.707` and the certificate would fail at `1.25`; the
per-level decay is what closes it. See [UFULL_ASSEMBLY_PROOF.md](UFULL_ASSEMBLY_PROOF.md).

Lemma C is now **proved** ([LEMMA_C_PROOF.md](LEMMA_C_PROOF.md)). The covector `c` is a 3x+1
exponential sum; via Lemma B's 2-adic shell method, the assembly-strength bound
`g_b < sqrt(3/4) = 0.866 < 0.961` follows, with the cross-scale identity `S_k(2 xi) = S_{k-1}(xi)` -
once the open link - now **derived** (Lemma H, an elementary parity split: one parity gives a geometric
series that vanishes, the other is exactly the next-scale sum). Lemma C inherits only the Half-Shift
Invariance dependence that Lemmas A and B already carry, and adds no new conditional input.

### What remains

The spectral-gap mathematics is complete and all-`k` (Coset-Uniformity proved in
[HALFSHIFT_S4_LEMMA_A_PROOF.md](HALFSHIFT_S4_LEMMA_A_PROOF.md); `cert(k) < 0.9005 < 1`, uniform). What
does **not** follow, and is **withdrawn**, is cycle elimination: the gap `=> no cycles` inference is
false ([CYCLE_CLAIM_REFUTED.md](CYCLE_CLAIM_REFUTED.md)). Ruling out non-trivial Collatz cycles would
require an instrument sensitive to the deterministic orbit (e.g. linear-forms-in-logs / height bounds,
as in classical cycle-length results), not this transfer-operator spectral gap. The only thing the
spectral-gap result still admits as future work is a Lean formalisation of the gap itself, which would
formalise a correct-but-not-cycle-eliminating statement.

---

## Reproduce everything

```bash
pip install -r requirements.txt
python generate_figures.py          # the four README figures, from the real operator
python audit_halfshift_s4.py        # coset-uniformity, S4, isometry parity-split + boundaries (0 violations)
python attack1_lemmaA_proof.py      # closed-form B*B = 2^-d I, all (a,b), to k=14
python lemmaB_fact1_rigorous.py     # collision bound coll <= 3*2^k, both parities, vs the true Syracuse fibre
python adv_tril_sep_correct.py      # ||tril(Q_D)|| matches the dense operator (<1e-10); chain to k=24
python verify_assembly.py           # the assembly: cert < 0.9005 < 1 (Lemma A+C); matches build_T to 6 digits
python explore_vb_profile.py        # the v_b profile and per-row S_a (Lemma C ground truth)
python probe_periodization.py       # g_b^2 = 2^b(2 h_b - h_{b+1}); sup = 9/16 (the Lemma C reduction)
python probe_gb_collision.py        # g_b^2 = E_p/4^p (single index, p=k-b); sum|S|^2 = 2^p coll(p)
python probe_autocorr.py            # E_p = 2^{p-1}(coll-A); Lemma C <=> coll(p)-A_p <= (9/8)2^p
python probe_shell_halfshift.py     # shell proof: j=0 cancels, diag=2^{p-1}, coll-A <= (3/2)2^p
python probe_homometry_proof.py     # Lemma H: parity split proving S_k(2 xi) = S_{k-1}(xi)
python verify_lemma_h.py            # independent integer-exact check of Lemma H (0 failures)
python probe_cycle_link.py          # the 3x+1 vs 3x-1 control: same gap, but 3x-1 HAS cycles (retraction)
python probe_cycle_recovery.py      # cycle-detector tests: spectrum/traces are cycle-blind
```

## Repository map

| File | Role |
|------|------|
| `HALFSHIFT_S4_LEMMA_A_PROOF.md` | Lemma A: within-level isometry + Half-Shift coset-uniformity (proof) |
| `LEMMA_B_PROOF.md` | Lemma B: the collision bound `coll <= 3*2^k` (proof) |
| `STEP4_BLOCK_FORMULA_FOUNDATION.md` | the operator split `U = U_clean + D` and block formula |
| `HalfShiftInvariance_DRAFT.md` | Half-Shift Invariance / rank-1 `S_odd` (Lean draft) |
| `CYCLE_CLAIM_REFUTED.md` | **the retraction**: why the spectral gap does not eliminate cycles (3x-1 control) |
| `CYCLE_STRUCTURE_RECOVERY.md` | follow-up: no cycle structure is recoverable from the spectral side; where it lives |
| `probe_cycle_link.py`, `probe_cycle_recovery.py` | the 3x+1 vs 3x-1 control + cycle-detector tests |
| `UFULL_ASSEMBLY_PROOF.md` | the row-sum assembly: `cert(k) < 0.9005 < 1` from Lemma A + Lemma C |
| `LEMMA_C_PROOF.md` | Lemma C proved (shell method + Lemma H homometry); assembly-strength bound |
| `UFULL_ASSEMBLY_PLAN.md` | the prior cold-start brief for the assembly (now actioned) |
| `verify_assembly.py`, `explore_vb_profile.py`, `verify_lemma_h.py`, `probe_*.py` | assembly + Lemma C / H verification |
| `analytic_proofs.py` | `build_T(k)`: the transfer operator (the numerical oracle) |
| `audit_halfshift_s4.py`, `attack1_lemmaA_proof.py`, `attack1_Sodd.py` | Lemma A / CU / S4 verification |
| `lemmaB_fact1_rigorous.py`, `attack3_*.py` | Lemma B verification |
| `adv_tril_sep_correct.py`, `step4_*.py`, `diagnose_block_formula.py` | block / defect numerics |
| `generate_figures.py` | regenerates the figures |

## References

- L. Mori (2024), *C\*(T_1,T_2) irreducible on l^2(N) iff Collatz*, arXiv:2411.08084.
- T. Tao (2019), *Almost all Collatz orbits attain almost bounded values*, arXiv:1909.03562.
- A. Kontorovich, J. Lagarias (2009-2010), transfer operators for 3x+1.

## License

Public domain. The detailed exploratory history (many dead ends) lives in a separate private
repository; this repo is the curated, reproducible result.

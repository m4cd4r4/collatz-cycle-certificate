# The spectral gap does NOT eliminate cycles (refutation of the central claim)

Written 2026-06-02 (Session 40). This document records a **negative result that retracts the
headline claim** of this repository. It was found while trying to write out the spectral reduction
`cert(k) < 1 => |lambda_2(U_k)| < 1 => no non-trivial cycles` in full. The first arrow is sound; the
**second arrow is false**, and a clean control experiment proves it.

## What the certificate actually establishes (still true)

Lemmas A, B, C and the assembly rigorously prove, for the Syracuse transfer operator `U_k` (the
column-stochastic mod-`2^k` Markov operator `T_k` in the character basis):
```
    cert(k) = max_a sum_b Q_k[a,b] 2^{a-b} < G_up + 2^{-3/2} = 0.9005 < 1,   uniformly in k,
```
hence `|lambda_2(U_k)| <= rho(Q_k) <= cert(k) < 1`: a **uniform spectral gap** for `T_k`. That part of
the program is correct and unaffected by what follows.

## The false inference

The repository claimed: a uniform spectral gap `|lambda_2(U_k)| < 1` **eliminates non-trivial
Syracuse cycles**. It does not. A transfer/Markov operator's spectral gap measures mixing of the
*lift-averaged* dynamics (each residue `r` is spread over all targets `Syr(r + m 2^k) mod 2^k`,
`m < 2^k`). Averaging over the lift destroys the deterministic orbit structure, so a genuine cycle
of the integer map leaves **no trace** in the spectrum of `T_k`. This is the same phenomenon as the
doubling map `x -> 2x mod 1`, which is exponentially mixing (spectral gap) yet has a dense set of
periodic points.

## The control experiment (decisive)

The `3x-1` map (odd core `n -> oddpart(3n-1)`) is identical to Syracuse except for the sign, so the
transfer operator is built the same way. **Unlike `3x+1`, it has known non-trivial cycles:**
`{5, 7}` (`5 -> 7 -> 5`) and `{17, 25, 37, 55, 41, ...}`. Building its operator and certificate with
the same code (`probe_cycle_link.py`, sign flip only):

| k | 3x+1 cert | 3x+1 rho(Q) | 3x+1 \|l2\| | 3x-1 cert | 3x-1 rho(Q) | 3x-1 \|l2\| |
|---|---|---|---|---|---|---|
| 6  | 0.6338 | 0.5535 | 0.277 | 0.6018 | 0.5039 | 0.297 |
| 8  | 0.6347 | 0.5661 | 0.255 | 0.6061 | 0.5272 | 0.293 |
| 10 | 0.6345 | 0.5673 | 0.270 | 0.6061 | 0.5365 | 0.291 |
| 11 | 0.6343 | 0.5664 | 0.252 | 0.6060 | 0.5388 | 0.283 |

The `3x-1` operator passes the certificate **more strongly** than `3x+1` (`cert ~ 0.606 < 0.634`,
`rho(Q) ~ 0.538 < 0.567`, `|l2| ~ 0.29 < 1`), uniformly in `k` - yet `3x-1` provably **has**
non-trivial cycles. Therefore:

> **`cert(k) < 1` (a uniform spectral gap of the transfer operator) does NOT imply the absence of
> non-trivial cycles.** The identical certificate machinery "certifies" a map that demonstrably has
> cycles. The inference the repository's conclusion rested on is invalid.

## What survives, and what is withdrawn

**Withdrawn:** the claim that this program proves "no non-trivial Syracuse/Collatz cycles." It does
not. The certificate is blind to cycles.

**Survives (as correct mathematics, now without the cycle corollary):**
- Lemma A: `||P_a U_clean P_b||_2 = 2^{-(b-a)/2}` exactly, all `k` (CU + SB + S4).
- Lemma B: `||tril(Q_k)||_2 <= sqrt(3) 2^{-k/2}`, all `k` (the collision bound).
- Lemma C: `v_b <= (3/4) 2^{-b} 2^{-k/2}`, all `k` (Lemma H homometry + shell bound).
- Coset-Uniformity (CU): proved, all `k`, elementary.
- The assembly: `cert(k) < 0.9005 < 1` uniformly.
These together are a correct, self-contained proof of a **uniform spectral gap for the Syracuse
transfer operator** `T_k`. That is a real (if modest, and likely known-in-spirit) result. It simply
does not have the cycle-elimination consequence that was claimed.

## Why the gap cannot see cycles (the mechanism, for the record)

A non-trivial cycle `{a_1, ..., a_L}` is a periodic orbit of the *deterministic* map. For the
operator `T_k` to "see" it, the cycle residues would have to form an invariant set under `T_k`. But
`T_k[s, r] = 2^{-k} #{m < 2^k : Syr(r + m 2^k) = s mod 2^k}` mixes each residue over a full coset
(coset-uniformity - the very lemma the gap proof relies on). The `m = 0` term carries the true edge
`a_i -> a_{i+1}` with weight only `2^{-k}`; the other `2^k - 1` lifts scatter elsewhere. So the cycle
is diluted to invisibility, and its existence is consistent with `|lambda_2| < 1`. A spectral gap of
the a.c./averaged transfer operator is simply the wrong instrument for periodic orbits.

## Honest status of the repository

The repository now stands as: a correct elementary proof of a uniform spectral gap for the Syracuse
mod-`2^k` transfer operator, **with the cycle-elimination claim retracted**. Ruling out non-trivial
Collatz cycles would require an instrument sensitive to the deterministic orbit structure (e.g. a
height/transcendence or linear-forms-in-logs argument as in the classical cycle bounds of Baker-type,
or a genuinely different operator), not the spectral gap proved here.

Reproduce: `python probe_cycle_link.py`.

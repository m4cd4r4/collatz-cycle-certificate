# -*- coding: utf-8 -*-
"""
CORRECTED separable ||tril(Q_D)||_2 to high k, plus the full bound chain.

Validated construction (matches dense build_T tril exactly at k<=12):
  u_a = ||P_a e_{r*}||  with u_a^2 = d_a/N   (d_a=2^{k-2-a}, N=2^{k-1})
  v_b = ||P_b c||,  c = D[r*,:] = counts/2^k  (state-basis covector),
        v_b^2 = sum_{xi: v2(xi)=b} |<chi_xi, c>|^2   computed via direct char projection.
  tril(Q_D)[a,b] = u_a v_b for a>=b.

The char projection <chi_xi, c> for xi=0..N-1 over odd residues n=2s+1 is:
  P[xi] = (1/sqrt(N)) sum_s c[s] w^{-(2s+1) xi},  w=exp(2pi i/2^k).
This is a length-N transform; do it with one FFT per k by noting
  sum_s c[s] w^{-(2s+1)xi} = w^{-xi} sum_s c[s] (w^{-2xi})^s.
For xi in 0..N-1, w^{-2xi} = exp(-2pi i (2xi)/2^k) = exp(-2pi i xi /2^{k-1}) = standard N-th root.
So sum_s c[s] (rootN)^{-xi s} = FFT_N(c)[xi]. Hence P[xi] = w^{-xi} FFT_N(c)[xi] / sqrt(N).
=> |P[xi]| = |FFT_N(c)[xi]| / sqrt(N).  Clean, O(N log N), memory O(N).

Then verify the bound chain:
  ||tril(Q_D)|| <= ||Q_D||_2 <= ||D||_2 = sqrt(coll)/2^k <= sqrt(3) 2^{-k/2}.
"""
import sys, os, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from analytic_proofs import build_T, v2, syr


def rstar(k):
    return (-pow(3, -1, 1 << k)) % (1 << k)


def covector_state(k):
    """c[s] = (1/2^k) #{m: Syr(r*+m2^k)=2s+1 mod 2^k},  s=0..N-1."""
    mod = 1 << k
    N = mod >> 1
    rs = rstar(k)
    c = np.zeros(N)
    for m in range(mod):
        t = syr(rs + m * mod) % mod
        c[(t - 1) // 2] += 1.0
    coll = int(np.sum(c.astype(np.int64) ** 2))
    return c / mod, coll


def level_energies(k):
    """v_b^2 = sum_{xi: v2(xi)=b, 1<=xi<N} |P[xi]|^2,  |P[xi]| = |FFT_N(c)[xi]|/sqrt(N)."""
    N = 1 << (k - 1)
    c, coll = covector_state(k)
    Pmag = np.abs(np.fft.fft(c)) / np.sqrt(N)     # length N, |<chi_xi,c>|
    e = np.zeros(k - 1)
    for xi in range(1, N):
        b = v2(xi)
        if b <= k - 2:
            e[b] += Pmag[xi] ** 2
    return e, coll, Pmag


def tril_QD(k):
    e, coll, Pmag = level_energies(k)
    N = 1 << (k - 1)
    u = np.array([np.sqrt((1 << (k - 2 - a)) / N) for a in range(k - 1)])
    v = np.sqrt(e)
    nlev = k - 1
    QD = np.zeros((nlev, nlev))
    for a in range(nlev):
        for b in range(nlev):
            if a >= b:
                QD[a, b] = u[a] * v[b]
    triln = np.linalg.svd(QD, compute_uv=False)[0]
    # full Q_D (all a,b) for ||Q_D||
    QDfull = np.outer(u, v)
    QDnorm = np.linalg.svd(QDfull, compute_uv=False)[0]
    return triln, QDnorm, coll


def tril_full_buildT(k):
    U = build_T(k).T
    mod = 1 << k
    N = mod >> 1
    odds = (2 * np.arange(N) + 1).astype(np.int64)
    w = np.exp(2j * np.pi / mod)
    C = (w ** (np.outer(odds, np.arange(N)))) / np.sqrt(N)
    Uc = C.conj().T @ U @ C
    lev = np.array([(v2(int(x)) if x != 0 else -1) for x in range(N)])
    levels = sorted(set(lev[lev >= 0].tolist()))
    idx = {a: np.where(lev == a)[0] for a in levels}
    Q = np.zeros((len(levels), len(levels)))
    for ia, a in enumerate(levels):
        for ib, b in enumerate(levels):
            blk = Uc[np.ix_(idx[a], idx[b])]
            if blk.size:
                Q[ia, ib] = np.linalg.svd(blk, compute_uv=False)[0]
    return np.linalg.svd(np.tril(Q), compute_uv=False)[0]


if __name__ == "__main__":
    sqrt3 = np.sqrt(3.0)
    print("CROSS-CHECK corrected separable tril(Q_D) vs dense build_T tril(Q_full), k<=12:")
    for k in range(6, 13):
        tD, _, _ = tril_QD(k)
        tF = tril_full_buildT(k)
        print(f"  k={k:2d}: tril_QD={tD:.8e}  tril_full={tF:.8e}  reldiff={abs(tD-tF)/tF:.2e}")
    print()
    print("HIGH-k bound chain:  ||tril(Q_D)|| <= ||Q_D|| <= ||D||=sqrt(coll)/2^k <= sqrt(3)2^-k/2")
    print(f"{'k':>3} {'tril_QD':>14} {'||Q_D||':>13} {'||D||=sqcoll/2^k':>16} {'sqrt3*2^-k/2':>13} "
          f"{'chain_ok':>9} {'tril*2^k/2':>11} {'||D||*2^k/2':>11}")
    all_ok = True
    for k in range(6, 25):
        tD, QDn, coll = tril_QD(k)
        Dn = np.sqrt(coll) / (1 << k)
        bound = sqrt3 * 2 ** (-k / 2)
        chain = (tD <= QDn * (1 + 1e-9)) and (QDn <= Dn * (1 + 1e-9)) and (Dn <= bound * (1 + 1e-9))
        all_ok &= chain
        print(f"{k:>3} {tD:>14.6e} {QDn:>13.6e} {Dn:>16.6e} {bound:>13.6e} "
              f"{str(chain):>9} {tD*2**(k/2):>11.6f} {Dn*2**(k/2):>11.6f}", flush=True)
    print(f"\n  ALL chain holds (tril<=||Q_D||<=||D||<=sqrt3*2^-k/2): {all_ok}")
    print(f"  sqrt(31/12)={np.sqrt(31/12):.6f} (limit of ||D||*2^k/2)  sqrt(3)={sqrt3:.6f}")

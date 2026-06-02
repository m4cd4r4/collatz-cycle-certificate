# -*- coding: utf-8 -*-
"""
DECISIVE TEST of the (L2) link: does |lambda_2(T_k)| < 1 imply "no non-trivial cycles"?

The 3x+1 (Syracuse) map has NO known non-trivial cycle. The 3x-1 map DOES have
non-trivial cycles: {5,7} (5 -> 7 -> 5) and {17,25,37,55,...}. If the analogous
transfer operator T_k for 3x-1 ALSO has |lambda_2| bounded below 1 uniformly,
then a uniform spectral gap does NOT preclude cycles, and the (L2) inference
"gap => no cycles" is FALSE as a general principle. If instead |lambda_2(3x-1)|
-> 1 (cycles show up in the spectrum), (L2) has a mechanism.

We build T_k exactly as build_T does, but with the +1 replaced by -1.
"""
import numpy as np


def syr_pm(n, sign):
    """Syracuse-type odd core: oddpart(3n+sign).  sign=+1 (Collatz) or -1 (3x-1)."""
    val = 3 * n + sign
    if val == 0:
        return 0
    while val % 2 == 0:
        val //= 2
    return val


def build_T_pm(k, sign):
    mod = 1 << k
    N = mod >> 1
    T = np.zeros((N, N))
    for src in range(N):
        r = 2 * src + 1
        for m in range(mod):
            n = r + m * mod
            s = syr_pm(n, sign) % mod
            tgt = (s - 1) // 2 if s % 2 == 1 else None
            if tgt is not None and 0 <= tgt < N:
                T[tgt, src] += 1.0
    T /= mod
    return T


def lambda2(T):
    ev = np.linalg.eigvals(T)
    mags = np.sort(np.abs(ev))[::-1]
    return mags[0], mags[1]   # top (Perron ~1) and second


if __name__ == "__main__":
    print("Known 3x-1 cycle check: 5 -> oddpart(3*5-1=14)=7 -> oddpart(3*7-1=20)=5. Cycle {5,7}.")
    print(f"  syr(5,-1)={syr_pm(5,-1)}, syr(7,-1)={syr_pm(7,-1)}  (should be 7,5)\n")
    print(f"{'k':>3} | {'+1 (Collatz)':>22} | {'-1 (has cycles {5,7}...)':>26}")
    print(f"{'':>3} | {'|l1|':>10} {'|l2|':>10} | {'|l1|':>11} {'|l2|':>13}")
    for k in range(4, 13):
        Tp = build_T_pm(k, +1)
        Tm = build_T_pm(k, -1)
        l1p, l2p = lambda2(Tp)
        l1m, l2m = lambda2(Tm)
        print(f"{k:>3} | {l1p:>10.6f} {l2p:>10.6f} | {l1m:>11.6f} {l2m:>13.6f}")
    print("\nIf |l2| for 3x-1 stays bounded < 1 like Collatz, then 'gap => no cycles' is FALSE:")
    print("a uniform spectral gap does NOT detect the known {5,7} (and larger) 3x-1 cycles.")

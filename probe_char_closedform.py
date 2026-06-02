# -*- coding: utf-8 -*-
"""
Hunt for a closed form of |<chi_xi, c>| = |FFT_N(c)[xi]|/sqrt(N).
Hypothesis families:
  (H1) depends only on v2(xi)         -> would make levels degenerate (FALSE, there's spread)
  (H2) ~ 2^{-v2(xi)} * (Gauss-sum factor depending on odd part)
  (H3) the LEVEL energy v_b^2 has a closed form even if per-char does not.

Key target: prove v_b^2 <= C^2 4^{-b} 2^{-k}  (Lemma C).
We test the cleaner equivalent:  E_b := v_b^2 * 2^k * 4^b  is bounded (== g_b^2).
And we test the cumulative tail  sum_{b'>=b} v_{b'}^2  vs  4^{-b} 2^{-k} (Lemma-B-localised).
"""
import numpy as np
from adv_tril_sep_correct import covector_state, level_energies
from analytic_proofs import v2


def coeff_table(k, show=True):
    N = 1 << (k - 1)
    c, coll = covector_state(k)
    F = np.fft.fft(c)
    Pmag = np.abs(F) / np.sqrt(N)
    # |<chi,c>|^2 * 2^k, indexed; look at structure vs xi odd part & v2
    if show:
        print(f"k={k} N={N} coll={coll}")
        # For each xi, factor xi = 2^b * m, m odd. tabulate scaled coeff.
        print(f"{'xi':>6} {'b=v2':>4} {'m=odd':>6} {'|P|^2*2^k*4^b':>16}")
        # sample: all xi with small odd part m in {1,3,5,7} across levels
        for m in [1, 3, 5, 7]:
            for b in range(k - 1):
                xi = (1 << b) * m
                if 1 <= xi < N:
                    val = Pmag[xi] ** 2 * (1 << k) * (4 ** b)
                    print(f"{xi:>6} {b:>4} {m:>6} {val:>16.6f}")
            print()
    return c, Pmag, coll


def cum_tail(k):
    e, coll, _ = level_energies(k)   # e[b] = v_b^2
    nlev = len(e)
    # E_b = v_b^2 * 2^k * 4^b  (== g_b^2);  and tail sums
    Eb = e * (1 << k) * (4.0 ** np.arange(nlev))
    return Eb, e, coll


if __name__ == "__main__":
    np.set_printoptions(linewidth=200, suppress=True, precision=5)
    coeff_table(12)
    print("=" * 70)
    print("E_b = v_b^2 * 2^k * 4^b  (= g_b^2), across k -- is sup_b E_b bounded?")
    print(f"{'k':>3}  {'max_b E_b':>10}  {'argmax b':>9}  {'E_0':>8}  {'E_b profile (sqrt=g_b)':>10}")
    for k in range(6, 23, 2):
        Eb, e, coll = cum_tail(k)
        print(f"{k:>3}  {Eb.max():>10.6f}  {int(Eb.argmax()):>9}  {Eb[0]:>8.6f}   g_b={np.sqrt(Eb)}")
    print("\nIf sup_b g_b is bounded by a universal C across all k, Lemma C holds with that C.")
    print("Candidate closed forms for limit g_b (bulk):  sqrt? ", np.sqrt(0.661807),
          "  0.6618^2=", 0.661807**2)

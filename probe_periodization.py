# -*- coding: utf-8 -*-
"""
Periodization handle on Lemma C.  Parseval/aliasing identity:
  sum_{xi ≡ 0 mod 2^b} |F[xi]|^2  =  (N/2^b) * sum_{r=0}^{N/2^b-1} |P_b c[r]|^2
where P_b c[r] = sum over the (2^b) aliases  c[r + t*(N/2^b)], t=0..2^b-1
(the 2^b-fold periodization of c), and F = FFT_N(c).

Hence the band energy at EXACTLY level b is
  B_b := sum_{xi: v2(xi)=b} |F[xi]|^2 = M2b - M2b1
where M_j := sum_{xi ≡ 0 mod 2^j}|F[xi]|^2 = (N/2^j) ||fold_{2^j}(c)||^2.
And  g_b^2 = 2 * 4^b * B_b   (since v_b^2 = B_b/N, g_b^2 = 2^k 4^b v_b^2 = 2*4^b B_b).

Goal: see that fold_{2^j}(c) has a clean closed form (it should: c is the r* fiber
distribution; its periodizations collapse to simple dyadic vectors), which would
prove g_b^2 <= 9/16.
"""
import numpy as np
from adv_tril_sep_correct import covector_state
from analytic_proofs import v2


def fold(c, twoj):
    """2^j-fold periodization: sum c[r + t*(N/2^j)] over t, length N/2^j."""
    N = len(c)
    M = N // twoj
    return c.reshape(twoj, M).sum(axis=0)


def band_energies(k):
    N = 1 << (k - 1)
    c, coll = covector_state(k)
    F = np.fft.fft(c)
    # direct band energies
    B = np.zeros(k - 1)
    for xi in range(1, N):
        b = v2(xi)
        if b <= k - 2:
            B[b] += abs(F[xi]) ** 2
    # via periodization: M_j = (N/2^j)||fold_{2^j}(c)||^2  (this counts xi≡0 mod 2^j incl xi=0)
    M = []
    for j in range(0, k):
        twoj = 1 << j
        if twoj > N:
            break
        fc = fold(c, twoj)
        M.append((N / twoj) * np.sum(fc ** 2))
    M = np.array(M)
    g2 = 2.0 * (4.0 ** np.arange(k - 1)) * B
    return B, M, g2, c


if __name__ == "__main__":
    np.set_printoptions(linewidth=200, suppress=True, precision=6)
    for k in [8, 10, 12]:
        B, M, g2, c = band_energies(k)
        print(f"\n=== k={k} ===")
        print(f"  g_b^2 (=2*4^b*B_b): {g2}")
        print(f"  max g_b^2 = {g2.max():.6f}  (9/16={9/16:.6f})")
        # B_b = M_b - M_{b+1}  exactly, for ALL b including b=0.  F[0]=1 (the xi=0 term) sits in
        # every M_j (0 == 0 mod 2^j for all j), so it cancels in every difference -- no correction.
        bandM = M[:-1] - M[1:]
        print(f"  B_b direct      : {B}")
        print(f"  M_b - M_{{b+1}}    : {bandM[:len(B)]}   (matches B_b?)")
        # the fold vectors: do they collapse to simple form?
    print("\nFold structure of c at k=12 (periodizations collapse?):")
    B, M, g2, c = band_energies(12)
    N = len(c)
    for j in range(0, 8):
        fc = fold(c, 1 << j)
        nz = np.count_nonzero(fc > 1e-12)
        print(f"  j={j}: len={len(fc):5d} nnz={nz:5d} sum={fc.sum():.6f} "
              f"max={fc.max():.6f} ||.||^2={np.sum(fc**2):.8e} distinctvals={len(set(np.round(fc,9)))}")

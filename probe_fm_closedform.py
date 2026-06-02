# -*- coding: utf-8 -*-
"""
Pin the odd-part Gauss factor f(m) and the exact boundary g_b values.
  Bulk hypothesis:  |<chi_{2^b m}, c>|^2 * 2^k * 4^b  ->  f(m)   (b away from top)
  Then  g_b^2 = sum_{m odd, 2^b m < N} f(m) * (boundary corrections near top).

Targets to explain exactly:
  bulk g_b^2 -> 7/16 = 0.4375          (so sum_{m odd>=1} f(m) = 7/16 ?)
  sup_b g_b^2 = 9/16  at b=k-4         (the sharp Lemma-C constant 3/4)
  tail g_{k-2}^2=1/2, g_{k-3}^2=1/8, g_{k-4}^2=9/16, g_{k-5}^2=0.3437..
"""
import numpy as np
from adv_tril_sep_correct import covector_state
from analytic_proofs import v2


def f_of_m(k, mmax=33):
    N = 1 << (k - 1)
    c, coll = covector_state(k)
    Pmag = np.abs(np.fft.fft(c)) / np.sqrt(N)
    out = {}
    for m in range(1, mmax, 2):
        # use a mid-level b to stay in the bulk (away from top boundary)
        b = max(0, (k - 1) // 3)
        xi = (1 << b) * m
        if xi < N:
            out[m] = Pmag[xi] ** 2 * (1 << k) * (4 ** b)
    return out, coll


def exact_boundary_gb(k):
    """g_b^2 for the top several levels, as exact-ish fractions."""
    from adv_tril_sep_correct import level_energies
    e, coll, _ = level_energies(k)
    nlev = len(e)
    Eb = e * (1 << k) * (4.0 ** np.arange(nlev))   # = g_b^2
    return Eb


if __name__ == "__main__":
    np.set_printoptions(linewidth=200, suppress=True, precision=6)
    k = 16
    fm, coll = f_of_m(k)
    print(f"k={k}  odd-part Gauss factor f(m) = |<chi_{{2^b m}},c>|^2 * 2^k * 4^b (bulk b):")
    tot = 0.0
    for m, val in fm.items():
        tot += val
        # try to recognise f(m): is it ~ 1/(something)?
        print(f"  m={m:>3}  f={val:.8f}   1/f={1/val if val>1e-12 else 0:.3f}   "
              f"f*2^?={val* (1<<4):.5f}   cumsum={tot:.6f}")
    print(f"  sum_m f(m) (m<{33}) = {tot:.6f}   target bulk g^2 = 7/16 = {7/16:.6f}")

    print("\nTop-level g_b^2 (boundary), should be k-independent fixed tail:")
    for k in [10, 12, 14, 16, 18]:
        Eb = exact_boundary_gb(k)
        tail = Eb[-6:]
        print(f"  k={k}: g_b^2 top6 = {tail}   (b=k-6..k-2)")
    print(f"\n  Recognise tail as fractions:")
    Eb = exact_boundary_gb(16)
    from fractions import Fraction
    for j, val in enumerate(Eb[-6:]):
        fr = Fraction(val).limit_denominator(256)
        print(f"    b=k-{6-j} (top idx {len(Eb)-6+j}): g^2={val:.6f} ~ {fr}")
    print(f"\n  9/16={9/16:.6f} (sup, at b=k-4)   7/16={7/16:.6f} (bulk)   "
          f"1/2={1/2:.4f} (b=k-2)   1/8={1/8:.4f} (b=k-3)")

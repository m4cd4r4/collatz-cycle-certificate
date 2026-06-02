# -*- coding: utf-8 -*-
"""
Pin the v_b profile and the lower+diagonal weighted row-sum S_a, exactly,
across k. Groundwork for the U_full assembly (UFULL_ASSEMBLY_PLAN.md step 1).

  v_b   = ||P_b c||           (level energies of the defect covector)
  u_a   = 2^{-(a+1)/2}        (since u_a^2 = d_a/N = 2^{-1-a})
  S_a   = sum_{b<=a} Q_D[a,b] 2^{a-b}
        = u_a 2^a sum_{b<=a} v_b 2^{-b}
        = 2^{(a-1)/2} * sum_{b<=a} v_b 2^{-b}
  certificate lower+diag row-sum = max_a S_a   (want < ~0.453, ideally -> 0)

Also reports the scaled profiles to find the closed form of v_b:
  v_b * 2^{k/2}            (does the profile collapse across k?)
  v_b * 2^{b}             (geometric-in-b hypothesis)
"""
import numpy as np
from adv_tril_sep_correct import level_energies


def profiles(k):
    e, coll, _ = level_energies(k)
    v = np.sqrt(e)                       # v_b, b = 0..k-2
    nlev = k - 1
    a_idx = np.arange(nlev)
    u = 2.0 ** (-(a_idx + 1) / 2.0)      # u_a
    # S_a = 2^{(a-1)/2} * cumsum_{b<=a} v_b 2^{-b}
    partial = np.cumsum(v * 2.0 ** (-a_idx))
    S = 2.0 ** ((a_idx - 1) / 2.0) * partial
    return v, u, S, coll


if __name__ == "__main__":
    np.set_printoptions(linewidth=200, suppress=True, precision=6)
    for k in [6, 8, 10, 12, 14, 16, 18, 20]:
        v, u, S, coll = profiles(k)
        b = np.arange(len(v))
        print(f"\n===== k={k}  coll={coll}  coll/2^k={coll/(1<<k):.5f} =====")
        print(f"  max_a S_a       = {S.max():.6f}   (argmax a={int(S.argmax())}, top row a={k-2})")
        print(f"  S_a at top row  = {S[-1]:.6f}")
        print(f"  ||c||^2=sum v^2 = {(v**2).sum():.6e}   (<= 3*2^-k = {3/(1<<k):.6e})")
        # profile shape diagnostics
        print(f"  v_b              : {v}")
        print(f"  v_b * 2^(k/2)    : {v * 2**(k/2)}")
        print(f"  v_b * 2^b        : {v * 2.0**b}")
        print(f"  v_b*2^b*2^(k/2)  : {v * 2.0**b * 2**(k/2)}")
        print(f"  S_a              : {S}")

# -*- coding: utf-8 -*-
"""
PROVABLE all-k UPPER BOUND on coll(k), sufficient for Lemma B (we need an upper bound, not the exact 31/12).

Setup (proven structure):
  cf[t] = #{m in Z/2^k : oddpart(1+3m) mod 2^k = t}.  coll = sum_t cf[t]^2.
  Valuation shells A_j = {x=1+3m : v2(x)=j}, |A_j| = 2^{k-1-j} for j=0..k-1, plus a single top atom
  (the unique x=1+3m with v2(x)>=k, which exists since 3 is a unit so 1+3m≡0 mod 2^k has one solution m).
  R_j = { (x/2^j) mod 2^k : x in A_j }  (the residues hit by shell j).
  FACT 1 (proven): within a shell the map x -> x/2^j mod 2^k is INJECTIVE, so |R_j| = |A_j| = 2^{k-1-j}
     and cf restricted to shell j is 0/1.  [reason: x/2^j runs over an AP of odds with step 3 (a unit)
     of length 2^{k-1-j} <= 2^{k-1}; an AP of distinct-step length <= half-period hits distinct residues.]
  Therefore cf[t] = #{ j : t in R_j } = number of shells whose residue-set contains t.

  coll = sum_t cf[t]^2 = sum_t (#shells hitting t)^2 = sum_{j,j'} |R_j ∩ R_{j'}|
       = sum_j |R_j|  +  2 sum_{j<j'} |R_j ∩ R_{j'}|.
  diag = sum_j |R_j| = (2^k - 1) + 1(top atom) = 2^k   [EXACT, PROVEN].
  cross <= sum_{j<j'} |R_{j'}| = sum_{j'} j' * |R_{j'}|  (j' has j' lower partners j=0..j'-1).
       sum_{j'=0}^{k-1} j' 2^{k-1-j'} <= 2^{k-1} sum_{j'>=0} j' 2^{-j'} = 2^{k-1} * 2 = 2^k.
  => coll <= 2^k + 2 * 2^k = 3 * 2^k.   PROVABLE for all k.   (measured limit 31/12 ~ 2.583 < 3.)

So  ||c^op||_2 = sqrt(coll)/2^k <= sqrt(3) 2^{-k/2}  PROVABLY for all k.
This is the clean rigorous Lemma-B scale bound.  (The sharp 31/12 is verified; 3 is provable.)

This script verifies FACT 1 (injectivity per shell) and the bound coll <= 3*2^k for all tested k.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from fractions import Fraction
from analytic_proofs import v2


def oddpart(x):
    while x % 2 == 0:
        x //= 2
    return x


if __name__ == "__main__":
    print("Verify FACT1 (per-shell injectivity), diag=2^k exact, and coll <= 3*2^k provable bound:")
    print(f"{'k':>3} {'coll':>10} {'coll/2^k':>9} {'diag':>8} {'diag=2^k?':>9} {'shell-inj?':>10} "
          f"{'cross<=sumj*|Rj|?':>16} {'<=3*2^k?':>9}")
    for k in range(3, 17):
        mod = 1 << k
        R = {}
        Acount = {}
        shell_injective = True
        for m in range(mod):
            x = 1 + 3 * m
            j = v2(x)
            Acount[j] = Acount.get(j, 0) + 1
            R.setdefault(j, []).append(oddpart(x) % mod)
        # injectivity: |set(R_j)| == |R_j list|
        for j, lst in R.items():
            if len(set(lst)) != len(lst):
                shell_injective = False
        Rs = {j: set(lst) for j, lst in R.items()}
        js = sorted(Rs)
        diag = sum(len(Rs[j]) for j in js)
        # exact coll
        from collections import Counter
        cf = Counter()
        for j in js:
            for t in Rs[j]:
                cf[t] += 1
        coll = sum(v * v for v in cf.values())
        # cross
        cross = 0
        for ia in range(len(js)):
            for ib in range(ia + 1, len(js)):
                cross += len(Rs[js[ia]] & Rs[js[ib]])
        cross_bound = sum(js.index(jp) * len(Rs[jp]) for jp in js)  # sum j'*|R_j'| with j'=rank
        # use actual j value not rank for partners count: number of lower shells = count of j<j'
        cross_bound = 0
        for ib in range(len(js)):
            jp = js[ib]
            cross_bound += ib * len(Rs[jp])  # ib lower shells
        print(f"{k:>3} {coll:>10} {coll/mod:>9.5f} {diag:>8} {str(diag==mod):>9} "
              f"{str(shell_injective):>10} {str(cross<=cross_bound):>16} {str(coll<=3*mod):>9}")

    print()
    print("=" * 80)
    print("The provable chain:  coll = diag + 2*cross,  diag=2^k,  cross <= sum_{j'} j'|R_j'| <= 2^k.")
    print("  => coll <= 3*2^k  =>  ||c^op|| = sqrt(coll)/2^k <= sqrt(3) 2^{-k/2}  (ALL k, rigorous).")
    print("  Sharp (verified k<=18):  coll = (31*2^k +-c)/12  => ||c^op|| -> sqrt(31/12) 2^{-k/2} = 1.607.")
    print("=" * 80)
    print(f"  sqrt(3) = {np.sqrt(3):.6f}   sqrt(31/12) = {np.sqrt(31/12):.6f}")
    # the analytic cross bound series
    s = sum(Fraction(jp, 1 << jp) for jp in range(0, 100))
    print(f"  sum_{{j>=0}} j/2^j = {float(s)} (=2 exactly) => cross <= 2^{{k-1}}*2 = 2^k. QED bound.")

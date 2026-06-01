# -*- coding: utf-8 -*-
"""
PROVE the collision closed form coll(k) = (31*2^k +- c)/12 from the EXACT fiber structure.

The r* fiber (k even, Q0=1):  t_m = Syr(r*+m2^k) mod 2^k = oddpart(1+3m) mod 2^k, m in Z/2^k.
Since 3 is a unit mod 2^k, x:=1+3m is uniform over Z/2^k as m ranges over Z/2^k.  BUT the Syracuse
ODD PART oddpart(x) for x with valuation 0 (x odd) is x itself, which can be ANY odd residue in
[1,2^k) -- so far matches uniform.  The subtlety the earlier (wrong-limit) script missed: oddpart(x)
mod 2^k is taken AFTER dividing out all 2s, and for x EVEN the result x/2^v is odd and < 2^{k-v},
but we still reduce mod 2^k (no-op since < 2^k).  So the EVEN-uniform model SHOULD be exact... yet it
gave limit 3, not 31/12.  RESOLVE the discrepancy: it must be that x = 1+3m is NOT uniform over all of
Z/2^k in the way that matters, OR Syr divides differently.  Investigate by DIRECT comparison of the
two histograms cf_rstar vs cf_uniform-oddpart, find where they differ.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from analytic_proofs import v2, syr


def rstar(k):
    return (-pow(3, -1, 1 << k)) % (1 << k)


def cf_rstar(k):
    mod = 1 << k; rs = rstar(k)
    cf = np.zeros(mod, dtype=np.int64)
    for m in range(mod):
        cf[syr(rs + m * mod) % mod] += 1
    return cf


def oddpart(x):
    if x == 0:
        return 0
    while x % 2 == 0:
        x //= 2
    return x


def cf_uniform_oddpart(k):
    """x = 1+3m mod 2^k, oddpart(x) mod 2^k. (the naive model)"""
    mod = 1 << k
    cf = np.zeros(mod, dtype=np.int64)
    for m in range(mod):
        x = (1 + 3 * m) % mod
        cf[oddpart(x) % mod] += 1
    return cf


def cf_true_syr_noreduce(k):
    """x=1+3m (NOT reduced mod 2^k first), Syr = oddpart(3*(r*+m2^k)+1) = oddpart(2^k(1+3m))=oddpart(1+3m)
    but 1+3m here is the REAL integer (m in 0..2^k-1), oddpart taken on the real integer then mod 2^k."""
    mod = 1 << k
    cf = np.zeros(mod, dtype=np.int64)
    for m in range(mod):
        x = 1 + 3 * m               # REAL integer, not reduced
        cf[oddpart(x) % mod] += 1
    return cf


if __name__ == "__main__":
    print("Compare three histograms; find the correct model that reproduces coll(k)=(31*2^k+8)/12:")
    for k in range(4, 13, 2):  # even k
        a = cf_rstar(k)
        b = cf_uniform_oddpart(k)
        c = cf_true_syr_noreduce(k)
        ca = int(np.sum(a.astype(object) ** 2))
        cb = int(np.sum(b.astype(object) ** 2))
        cc = int(np.sum(c.astype(object) ** 2))
        print(f"  k={k:2d}: coll_rstar={ca}  coll_unifmod={cb}  coll_noreduce={cc}  "
              f"rstar==noreduce? {np.array_equal(a,c)}")

    print()
    print("So the correct model is x=1+3m as a REAL integer in [1, 1+3(2^k-1)], oddpart then mod 2^k.")
    print("x ranges over the AP {1,4,7,...} of length 2^k, i.e. {x : x≡1 mod 3, 1<=x<=3*2^k-2}.")
    print("Decompose coll by v2(x): for x≡1 mod3, the valuation v2(x) and the odd-part-mod-2^k collision.")
    print()
    # decompose collision by (v_x, v_x') and the value of oddpart mod 2^k
    for k in [6, 8, 10]:
        mod = 1 << k
        xs = [1 + 3 * m for m in range(mod)]
        # group by oddpart mod 2^k
        from collections import Counter
        cnt = Counter(oddpart(x) % mod for x in xs)
        coll = sum(v * v for v in cnt.values())
        # by valuation of x
        valcnt = Counter(v2(x) for x in xs)
        print(f"  k={k:2d}: coll={coll}  (=(31*2^k+8)/12={(31*mod+8)//12})  "
              f"val-distribution of x=1+3m: {dict(sorted(valcnt.items()))}")
        # The AP x=1+3m: how many have v2(x)=j? x even <=> 3m odd <=> m odd; etc. Should be ~2^k/2^{j+1}.
        print(f"        normalized val freq (x2^k): "
              f"{ {j: round(c/mod*mod) for j,c in sorted(valcnt.items())} }")

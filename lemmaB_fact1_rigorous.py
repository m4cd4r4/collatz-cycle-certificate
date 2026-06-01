# -*- coding: utf-8 -*-
"""
RIGOROUS proof of FACT 1 (per-shell injectivity) -- the one load-bearing gap in
attack3_provable_bound.py's all-k proof of  coll(k) <= 3*2^k  (hence Lemma B).

Model: the r* fiber value is  t_m = Syr(r* + m 2^k) mod 2^k = oddpart(a + 3m) mod 2^k,
m = 0..2^k-1, where 3r*+1 = 2^{2 ceil(k/2)} so a = 2^{2 ceil(k/2) - k} in {1, 2}:
    a = 1  for k EVEN  (3r*+1 = 2^k),   a = 2  for k ODD  (3r*+1 = 2^{k+1}).
x = a+3m ranges over the REAL-integer AP {a, a+3, ..., a+3(2^k-1)} (NOT reduced mod 2^k
before taking the odd part). cf[t] = #{m : t_m = t}, coll = sum_t cf[t]^2. The model is
cross-checked against the true Syracuse fiber (syr) for every k, both parities.

Shell A_j = {x = 1+3m : v2(x) = j}. For x in A_j, oddpart(x) = u = x/2^j (odd).
The within-shell image is  u mod 2^k.

CLAIM (FACT 1): x -> (x/2^j) mod 2^k is injective on A_j, for every j = 0..k-1.

PROOF (the mechanism this script verifies):
  Suppose x = 2^j u, x' = 2^j u' in A_j with u == u' (mod 2^k), x != x'.
  (1) u == u' (mod 2^k)  =>  u - u' = 2^k s  for some integer s != 0
                         =>  x - x' = 2^j (u-u') = 2^(k+j) s.            [identity I]
  (2) x == x' == a (mod 3)  =>  3 | (x - x') = 2^(k+j) s.  gcd(3,2)=1  =>  3 | s.  [identity II]
  (3) x, x' in [a, a + 3*2^k - 3] subset [1, 3*2^k)  =>  |x - x'| < 3*2^k  =>  2^(k+j)|s| < 3*2^k
                              =>  2^j |s| < 3.                            [range R]
  s != 0 and 3 | s  =>  |s| >= 3  =>  2^j |s| >= 3, contradicting (3).  Hence s = 0, x = x'. QED

This script verifies, to high k and for EVERY shell, the three load-bearing facts
(I), (II), (R) and the resulting injectivity, plus |A_j| = 2^(k-1-j), diag = 2^k,
cross <= 2^k, coll <= 3*2^k.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction
from collections import Counter, defaultdict
from analytic_proofs import v2, syr


def oddpart(x):
    while x % 2 == 0:
        x //= 2
    return x


def ap_constant(k):
    """a = 2^{2 ceil(k/2) - k} in {1,2}: 1 for even k, 2 for odd k."""
    return 1 if k % 2 == 0 else 2


def coll_syr(k):
    """Ground-truth collision count from the true Syracuse fiber Syr(r*+m 2^k)."""
    mod = 1 << k
    rs = (-pow(3, -1, mod)) % mod
    cf = Counter()
    for m in range(mod):
        cf[syr(rs + m * mod) % mod] += 1
    return sum(c * c for c in cf.values())


def check_k(k):
    mod = 1 << k
    a = ap_constant(k)
    shells = defaultdict(list)           # j -> list of (x, u=x/2^j)
    for m in range(mod):
        x = a + 3 * m                    # REAL integer, not reduced
        j = v2(x)
        shells[j].append((x, x >> j))

    report = {"k": k, "a": a}

    # |A_j| = 2^{k-1-j} for j < k, plus single top atom (v2 >= k)
    sizes_ok = True
    top_atoms = 0
    for j, lst in shells.items():
        if j >= k:
            top_atoms += len(lst)
        elif len(lst) != (1 << (k - 1 - j)):
            sizes_ok = False
    report["sizes_ok"] = sizes_ok
    report["top_atoms"] = top_atoms      # expect exactly 1

    # FACT 1: within-shell injectivity of u mod 2^k, and the proof mechanism (I),(II),(R)
    inj_ok = True
    identity_I_ok = True
    identity_II_ok = True
    range_R_ok = True
    for j, lst in shells.items():
        if j >= k:
            continue
        seen = {}
        for (x, u) in lst:
            r = u % mod
            if r in seen:
                inj_ok = False           # a genuine within-shell collision (must never happen)
            seen[r] = x
        # verify the mechanism on EVERY ordered within-shell pair: if images were equal
        # the witness s would be a nonzero multiple of 3 with 2^j|s|<3 (impossible).
        L = lst
        for ia in range(len(L)):
            xa, ua = L[ia]
            for ib in range(ia + 1, len(L)):
                xb, ub = L[ib]
                diff = xa - xb
                # (R): range bound holds for the AP
                if not (abs(diff) < 3 * mod):
                    range_R_ok = False
                # (II): both == a mod 3 so 3 | diff
                if (xa % 3 != a) or (xb % 3 != a) or (diff % 3 != 0):
                    identity_II_ok = False
                # (I): if u==u' mod 2^k then diff == 2^{k+j} s. Equivalent: 2^{k+j} | diff.
                #     Verify the contrapositive structure: diff is divisible by 2^j exactly j-shell,
                #     and IF 2^k | (ua-ub) then 2^{k+j} | diff. Check the implication directly.
                if (ua - ub) % mod == 0:
                    if diff % (1 << (k + j)) != 0:
                        identity_I_ok = False
        # keep the inner double loop affordable: only run it fully for small shells
    report["inj_ok"] = inj_ok
    report["identity_I_ok"] = identity_I_ok
    report["identity_II_ok"] = identity_II_ok
    report["range_R_ok"] = range_R_ok

    # diag = sum_j |R_j| = 2^k exact ; cross <= 2^k ; coll <= 3*2^k
    Rsets = {j: set(u % mod for (_, u) in lst) for j, lst in shells.items()}
    diag = sum(len(s) for s in Rsets.values())
    cf = Counter()
    for s in Rsets.values():
        for t in s:
            cf[t] += 1
    coll = sum(c * c for c in cf.values())
    cross = (coll - diag) // 2
    # analytic cross bound  sum_{j'} (#lower shells)*|R_j'|
    js = sorted(Rsets)
    cross_bound = sum(js.index(jp) * len(Rsets[jp]) for jp in js)
    # ground truth: the shell model must reproduce the real Syracuse fiber's coll
    coll_true = coll_syr(k)
    report.update(diag=diag, coll=coll, cross=cross, cross_bound=cross_bound,
                  diag_is_2k=(diag == mod), cross_le_bound=(cross <= cross_bound),
                  cross_bound_le_2k=(cross_bound <= mod), coll_le_3_2k=(coll <= 3 * mod),
                  coll_true=coll_true, model_matches_syr=(coll == coll_true))
    return report


def frobenius_chain(k, coll):
    """Verify  ||u||^2 = 1 - 2^{1-k} < 1  and  ||tril(Q_D)||_F <= ||u||*||c|| <= sqrt3*2^{-k/2}."""
    N = 1 << (k - 1)
    # ||u||^2 = sum_{a=0}^{k-2} 2^{k-2-a}/N
    u2 = sum(Fraction(1 << (k - 2 - a), N) for a in range(k - 1))
    c_norm2 = Fraction(coll, 1 << (2 * k))           # ||c||^2 = coll/4^k
    frob_bound2 = u2 * c_norm2                        # (||u||*||c||)^2  >= ||tril||_F^2
    sqrt3_bound2 = Fraction(3, 1 << k)                # (sqrt3 * 2^{-k/2})^2 = 3/2^k
    return u2, c_norm2, frob_bound2, sqrt3_bound2


if __name__ == "__main__":
    print("FACT 1 rigorous-proof verification (mechanism I/II/R) + chain to coll<=3*2^k:\n")
    KMAX_FULL = 13   # full O(|A_j|^2) pair check up to here
    print(f"{'k':>3} {'a':>2} {'sizes':>6} {'top':>4} {'inj':>5} {'(I)':>5} {'(II)':>5} {'(R)':>5} "
          f"{'diag=2^k':>9} {'cross<=bd':>9} {'bd<=2^k':>8} {'coll':>9} {'=syr?':>6} {'coll<=3*2^k':>11}")
    for k in range(3, KMAX_FULL + 1):
        r = check_k(k)
        print(f"{k:>3} {r['a']:>2} {str(r['sizes_ok']):>6} {r['top_atoms']:>4} {str(r['inj_ok']):>5} "
              f"{str(r['identity_I_ok']):>5} {str(r['identity_II_ok']):>5} {str(r['range_R_ok']):>5} "
              f"{str(r['diag_is_2k']):>9} {str(r['cross_le_bound']):>9} {str(r['cross_bound_le_2k']):>8} "
              f"{r['coll']:>9} {str(r['model_matches_syr']):>6} {str(r['coll_le_3_2k']):>11}")

    print("\nFrobenius chain  ||tril(Q_D)||_F <= ||u||*||c|| <= sqrt3*2^{-k/2}, and ||u||^2<1:")
    print(f"{'k':>3} {'||u||^2':>14} {'<1?':>4} {'(||u||||c||)^2':>16} {'(sqrt3 2^-k/2)^2':>17} {'chain?':>7}")
    for k in range(4, KMAX_FULL + 1):
        r = check_k(k)
        u2, c2, fb2, b2 = frobenius_chain(k, r["coll"])
        print(f"{k:>3} {float(u2):>14.10f} {str(u2 < 1):>4} {float(fb2):>16.8e} "
              f"{float(b2):>17.8e} {str(fb2 <= b2):>7}")
    print("\n  ||u||^2 = 1 - 2^{1-k}  (verify closed form):")
    for k in range(4, 11):
        u2, *_ = frobenius_chain(k, 1)[:1] + (None,)
        u2 = sum(Fraction(1 << (k - 2 - a), 1 << (k - 1)) for a in range(k - 1))
        print(f"    k={k:2d}: ||u||^2 = {u2} ; 1 - 2^(1-k) = {1 - Fraction(1, 1 << (k - 1))} ; "
              f"match={u2 == 1 - Fraction(1, 1 << (k - 1))}")

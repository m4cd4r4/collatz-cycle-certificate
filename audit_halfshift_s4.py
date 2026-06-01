# -*- coding: utf-8 -*-
"""
AUDIT of the two load-bearing inputs after the pivot:
  (1) coset-uniformity (S1 = Half-Shift step 3): as m ranges over Z/2^K, Syr(x+m*2^K) mod 2^K
      is uniform over the coset q0 + 2^{K-v} Z, multiplicity 2^{K-v}.  (K=k+1 in the draft.)
  (2) S4 dead-band: Sodd(alpha,m)=0 iff j<=v2(alpha)<=k-2 (m=k-j); full 2^{m-1} iff alpha in {0,2^{k-1}}.
  (3) only-j=d-survives + isometry B*B=2^-d I -- PARITY-SPLIT (Lemma B's bug was odd k) and BOUNDARY
      blocks (a=0,b=1), (a=0,b=k-2), (a=b-1).
Read-only. PYTHONIOENCODING=utf-8.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from collections import Counter


def v2(n):
    if n == 0:
        return 10**9
    c = 0
    while n % 2 == 0:
        n //= 2; c += 1
    return c


def syr(n):
    t = 3 * n + 1
    while t % 2 == 0:
        t //= 2
    return t


# ---------- (1) coset-uniformity, K = k+1 as in the draft, all odd x with v=v2(3x+1)<k ----------
def check_coset_uniformity(k):
    K = k + 1
    modK = 1 << K
    worst = []
    vmax_seen = -1
    for x in range(1, 1 << k, 2):          # odd x < 2^k
        v = v2(3 * x + 1)
        if v >= k:                          # the r* exception; skip (handled separately)
            continue
        vmax_seen = max(vmax_seen, v)
        q0 = (3 * x + 1) >> v
        hits = Counter(syr(x + m * modK) % modK for m in range(modK))
        # claim: support = {(q0 + t*2^{K-v}) % modK : t} , each multiplicity 2^{K-v}
        step = 1 << (K - v)
        expected_support = {(q0 + t * step) % modK for t in range(1 << v)}
        mult_ok = all(c == (1 << (K - v)) for c in hits.values())
        supp_ok = (set(hits) == expected_support)
        if not (mult_ok and supp_ok):
            worst.append((x, v, supp_ok, mult_ok))
    return vmax_seen, worst


# ---------- (2) S4 dead-band, with boundary m=2 emphasis ----------
def check_s4(k):
    mod = 1 << k
    w = np.exp(2j * np.pi / mod)
    def Sodd(alpha, m):
        return sum(w ** ((alpha * q) % mod) for q in range(1 << m) if q % 2 == 1)
    fails = []
    for m in range(2, k):                   # m = k-j, j=1..k-2  => m=2..k-1
        for alpha in range(mod):
            S = abs(Sodd(alpha, m))
            va = v2(alpha)
            in_band = (k - m) <= va <= (k - 2)
            sharp = (alpha == 0) or (alpha == (1 << (k - 1)))
            if in_band and S > 1e-7:
                fails.append(('band-nonzero', m, alpha, va, S))
            if sharp and abs(S - (1 << (m - 1))) > 1e-7:
                fails.append(('sharp-notfull', m, alpha, va, S))
            if (not in_band) and (not sharp) and S < 1e-7:
                fails.append(('unexpected-zero', m, alpha, va, S))
    return fails


# ---------- (3) isometry via closed form, parity-split + boundary blocks ----------
def lemmaA_isom(k):
    mod = 1 << k; N = mod >> 1; inv3 = pow(3, -1, mod); w = np.exp(2j * np.pi / mod)
    def Sodd_closed(alpha, m):
        L = 1 << (m - 1); a = alpha % mod
        if (2 * a) % mod == 0:
            return L * w ** a
        return w ** a * (w ** ((a * (1 << m)) % mod) - 1) / (w ** ((2 * a) % mod) - 1)
    def hat(eta, xi, b):
        return sum(w ** ((xi * inv3) % mod) * Sodd_closed((eta - xi * (1 << j) * inv3) % mod, k - j)
                   for j in range(1, min(b, k - 1) + 1))
    worst = 0.0; boundary = {}
    for b in range(1, k - 1):
        for a in range(0, b):
            d = b - a
            etas = [x for x in range(N) if v2(x) == b]
            xis = [x for x in range(N) if v2(x) == a]
            H = np.array([[hat(eta, xi, b) for xi in xis] for eta in etas])
            G = (H.conj() @ H.T) / (N * N)
            err = float(np.max(np.abs(G - 2.0 ** (-d) * np.eye(len(etas)))))
            worst = max(worst, err)
            if (a, b) in [(0, 1), (0, k - 2), (b - 1, b)]:
                boundary[(a, b)] = err
    return worst, boundary


if __name__ == "__main__":
    print("=== (1) COSET-UNIFORMITY (K=k+1), all odd x with v2(3x+1)<k ===")
    print(f"{'k':>3} {'parity':>6} {'max v seen':>10} {'violations':>11}")
    for k in range(5, 12):
        vmax, bad = check_coset_uniformity(k)
        print(f"{k:>3} {'even' if k%2==0 else 'odd':>6} {vmax:>10} {len(bad):>11}  {bad[:2]}")

    print("\n=== (2) S4 DEAD-BAND (all alpha, m=2..k-1; m=2 is the b=k-2 boundary) ===")
    print(f"{'k':>3} {'parity':>6} {'fails':>6}")
    for k in range(6, 13):
        fa = check_s4(k)
        print(f"{k:>3} {'even' if k%2==0 else 'odd':>6} {len(fa):>6}  {fa[:2]}")

    print("\n=== (3) LEMMA A ISOMETRY (closed form), PARITY SPLIT + boundary blocks ===")
    print(f"{'k':>3} {'parity':>6} {'worst||G-2^-d I||':>18}  boundary (a,b):err")
    for k in range(6, 15):
        worst, bd = lemmaA_isom(k)
        bstr = "  ".join(f"{ab}:{e:.1e}" for ab, e in bd.items())
        print(f"{k:>3} {'even' if k%2==0 else 'odd':>6} {worst:>18.3e}  {bstr}")

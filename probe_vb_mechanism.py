# -*- coding: utf-8 -*-
"""
Probe the per-level decay v_b ~ C 2^{-b} 2^{-k/2}.
Goal: find the closed form / mechanism for g_b := v_b 2^b 2^{k/2} being bounded.

v_b^2 = sum_{xi: v2(xi)=b, 1<=xi<N} |P[xi]|^2,   |P[xi]| = |FFT_N(c)[xi]| / sqrt(N).
At level b there are exactly (#xi in [1,N) with v2(xi)=b) characters.
We look at the per-character magnitudes |P[xi]| within each level to see if they
are ~equidistributed (=> v_b^2 ~ n_b * (typical |P|^2)) and find the closed form.
"""
import numpy as np
from adv_tril_sep_correct import covector_state, rstar
from analytic_proofs import v2


def per_level_detail(k):
    N = 1 << (k - 1)
    c, coll = covector_state(k)
    Pmag = np.abs(np.fft.fft(c)) / np.sqrt(N)
    # group by level
    levels = {}
    for xi in range(1, N):
        b = v2(xi)
        if b <= k - 2:
            levels.setdefault(b, []).append(Pmag[xi])
    rows = []
    for b in sorted(levels):
        arr = np.array(levels[b])
        n_b = len(arr)
        vb2 = (arr ** 2).sum()
        vb = np.sqrt(vb2)
        rows.append((b, n_b, vb, arr.max(), arr.min(), arr.mean(),
                     vb * 2.0 ** b * 2 ** (k / 2)))
    return rows, coll, c, Pmag


if __name__ == "__main__":
    np.set_printoptions(linewidth=200, suppress=True, precision=6)
    k = 12
    rows, coll, c, Pmag = per_level_detail(k)
    print(f"k={k}  coll={coll}  N={1<<(k-1)}")
    print(f"  c[0:8] = {c[:8]}   (covector head; c = D[r*,:])")
    print(f"  rstar = {rstar(k)}")
    print(f"\n{'b':>3} {'n_b=#chi':>9} {'v_b':>12} {'maxP':>12} {'minP':>12} {'meanP':>12} {'g_b=v_b 2^b 2^k/2':>18}")
    for b, n_b, vb, mx, mn, mean, g in rows:
        print(f"{b:>3} {n_b:>9} {vb:>12.6e} {mx:>12.6e} {mn:>12.6e} {mean:>12.6e} {g:>18.6f}")
    # n_b should be 2^{k-2-b}. typical |P| per char at level b?
    print("\n  Check: n_b vs 2^{k-2-b}, and per-char energy vb^2/n_b:")
    for b, n_b, vb, mx, mn, mean, g in rows:
        print(f"    b={b:2d}: n_b={n_b:6d}  2^(k-2-b)={1<<max(k-2-b,0):6d}  "
              f"vb^2/n_b={vb**2/n_b:.6e}  sqrt={np.sqrt(vb**2/n_b):.6e}  "
              f"*2^(k-1)={vb**2/n_b*(1<<(k-1)):.4f}")
    # What is the largest |P[xi]| overall and where?
    big = np.argsort(Pmag[1:])[::-1][:6] + 1
    print(f"\n  top |P[xi]| chars: xi={big.tolist()}  v2={[v2(int(x)) for x in big]}  |P|={Pmag[big]}")
    print(f"  |P[xi]| for xi=2^j (pure levels): ")
    for j in range(k - 1):
        xi = 1 << j
        if xi < (1 << (k - 1)):
            print(f"    xi=2^{j}={xi:6d}  v2={v2(xi)}  |P|={Pmag[xi]:.6e}  *sqrt(N)={Pmag[xi]*np.sqrt(1<<(k-1)):.6f}")

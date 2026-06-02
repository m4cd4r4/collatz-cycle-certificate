# -*- coding: utf-8 -*-
"""
Can any cycle-specific structure be recovered from the mod-2^k framework?
Test candidate cycle-sensitive invariants on 3x+1 (no cycles) vs 3x-1 (cycles {5,7},{17,...}).

A genuine cycle detector must DISTINGUISH the two maps in a way that (a) reflects cycles and
(b) scales (does not just brute-force search). We test:

  (0) ground truth: brute-force cycles of the deterministic map within [1, 2^k).
  (A) full spectrum: eigenvalues of T_k near roots of unity (period-L cycle -> e^{2pi i/L}?).
  (B) traces tr(T_k^n): closed-path weight; periodic-orbit signal?
  (C) the DETERMINISTIC partial map D_k on [1,2^k) (Syr if image < 2^k): its cycles ARE the
      real cycles. This is cycle-sensitive but is brute force in disguise; included as the
      "cycle-sensitive => intractable" baseline.
  (D) Perron eigenvalue multiplicity / second-eigenvalue gap: does multi-cycle 3x-1 differ?
"""
import numpy as np
from probe_cycle_link import syr_pm, build_T_pm


def brute_cycles(k, sign, maxstart=None):
    """Deterministic cycles among odd n in [1, 2^k): iterate Syr until repeat or leave range."""
    mod = 1 << k
    lim = maxstart or mod
    seen_cycles = set()
    cycles = []
    for start in range(1, lim, 2):
        n = start
        path = []
        steps = 0
        while n < mod and n >= 1 and steps < 10000:
            path.append(n)
            n = syr_pm(n, sign)
            steps += 1
            if n in path:
                # found a cycle
                i = path.index(n)
                cyc = tuple(sorted(path[i:]))
                if cyc not in seen_cycles and len(cyc) >= 1:
                    seen_cycles.add(cyc)
                    cycles.append(cyc)
                break
    # keep only genuine cycles (closed under Syr within range)
    real = []
    for cyc in cycles:
        if all(syr_pm(x, sign) in cyc for x in cyc):
            real.append(cyc)
    return sorted(real, key=len)


def spectrum_near_roots(T, Lmax=8):
    ev = np.linalg.eigvals(T)
    # distance of each eigenvalue to nearest L-th root of unity, L=2..Lmax, excluding lambda=1
    out = {}
    for L in range(2, Lmax + 1):
        roots = np.exp(2j * np.pi * np.arange(L) / L)
        # nearest non-trivial eigenvalue (|.|>0.5) to any primitive-ish root != 1
        best = 1e9
        for lam in ev:
            if abs(lam - 1) < 1e-6:
                continue
            d = min(abs(lam - z) for z in roots if abs(z - 1) > 1e-9)
            best = min(best, d)
        out[L] = best
    return out


if __name__ == "__main__":
    print("=== (0) GROUND TRUTH: deterministic cycles within [1,2^k) ===")
    for k in [7, 9, 11]:
        cp = brute_cycles(k, +1)
        cm = brute_cycles(k, -1)
        print(f"  k={k}: 3x+1 cycles={cp}")
        print(f"        3x-1 cycles={cm}")

    print("\n=== (A) eigenvalues near L-th roots of unity (period-L signal?) ===")
    for k in [8, 10]:
        Tp = build_T_pm(k, +1); Tm = build_T_pm(k, -1)
        sp = spectrum_near_roots(Tp); sm = spectrum_near_roots(Tm)
        print(f"  k={k} 3x+1 min-dist-to-root_L: {{L:round(d,3) for...}}")
        print(f"        " + "  ".join(f"L{L}:{sp[L]:.3f}" for L in sp))
        print(f"  k={k} 3x-1 (has 2-cycle {{5,7}} -> expect eigenvalue near -1?):")
        print(f"        " + "  ".join(f"L{L}:{sm[L]:.3f}" for L in sm))

    print("\n=== (B) traces tr(T_k^n), n=1..6 (closed-path weight) ===")
    for k in [8, 10]:
        Tp = build_T_pm(k, +1); Tm = build_T_pm(k, -1)
        tp = [np.trace(np.linalg.matrix_power(Tp, n)).real for n in range(1, 7)]
        tm = [np.trace(np.linalg.matrix_power(Tm, n)).real for n in range(1, 7)]
        print(f"  k={k} 3x+1 tr(T^n): {[round(x,4) for x in tp]}")
        print(f"  k={k} 3x-1 tr(T^n): {[round(x,4) for x in tm]}")

    print("\n=== (C) DETERMINISTIC partial map D_k cycle count (cycle-sensitive baseline) ===")
    for k in [9, 11, 13]:
        cp = brute_cycles(k, +1); cm = brute_cycles(k, -1)
        print(f"  k={k}: 3x+1 #cycles={len(cp)} (sizes {[len(c) for c in cp]}); "
              f"3x-1 #cycles={len(cm)} (sizes {[len(c) for c in cm]})")

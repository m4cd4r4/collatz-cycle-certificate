# -*- coding: utf-8 -*-
"""
WITHIN-LEVEL structure of the clean upper block  B = P_a U_clean P_b,  and the EXACT mechanism that
makes it a SCALED PARTIAL ISOMETRY with common singular value 2^{-(b-a)/2}.

Plan:
  (1) Show B^* B = 2^{-(b-a)} * Pi  for an orthogonal projection Pi (=> all nonzero sv equal 2^{-d/2}).
      Equivalently  B B^* = 2^{-(b-a)} * (orthogonal projection).  This is the Lemma-A core identity.
  (2) Decompose the lower-triangular back-flow:  ||P_a U_full P_b||_2 for a>=b is ENTIRELY the r*
      defect (P_a U_clean P_b = 0 for a>=b since U_clean is strictly upper in level), so
        tril(Q_full) = level-matrix of the rank-1 defect, ||.||_2 = beta 2^{-k/2}, beta -> ~0.684.
  (3) Energy/norm reconciliation closed arithmetic:
        Q_mass[a,b] = (1/d_b)||B||_F^2 = (1/d_b)*sum sv^2 = (1/d_b)*rank*2^{-d} = (d_a/d_b)*2^{-d}? no.
        rank(B) = d_b (full source) for a<b, and d_a/d_b = 2^{?}. Pin it exactly.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from analytic_proofs import build_T, v2

np.set_printoptions(suppress=True, precision=6, linewidth=200)


def odds(k):
    return np.array([2 * s + 1 for s in range(1 << (k - 1))], dtype=np.int64)


def char_matrix(k):
    mod = 1 << k; N = mod >> 1
    w = np.exp(2j * np.pi / mod)
    return (w ** (np.outer(odds(k), np.arange(N)))) / np.sqrt(N)


def level_of_xi(k):
    N = 1 << (k - 1)
    return np.array([(v2(int(x)) if x != 0 else -1) for x in range(N)])


def U_clean_matrix(k):
    o = odds(k)
    v = np.array([v2(3 * int(r) + 1) for r in o])
    Uc = build_T(k).T.copy()
    Uc[v >= k, :] = 0.0
    return Uc


def bases_dict(k):
    C = char_matrix(k); lev = level_of_xi(k)
    mz = sorted(set(int(l) for l in lev if l >= 0))
    return {m: C[:, lev == m] for m in mz}, mz


if __name__ == "__main__":
    print("=" * 100)
    print("(1) LEMMA A CORE IDENTITY:  B = P_a U_clean P_b (a<b).  Is  B B^* = 2^{-(b-a)} * projection?")
    print("    (all nonzero sv = 2^{-d/2}  <=>  B B^* is 2^{-d} times an orthogonal projection)")
    print("=" * 100)
    for k in [8, 10, 12]:
        Uc = U_clean_matrix(k)
        bs, mz = bases_dict(k)
        print(f"\n  k={k}:")
        for (a, b) in [(0, 1), (0, 2), (1, 3), (2, 4), (0, 4)]:
            if a in bs and b in bs and a < b:
                Va, Vb = bs[a], bs[b]
                B = Va.conj().T @ (Uc @ Vb)         # d_a x d_b
                d = b - a
                G = B @ B.conj().T                  # d_a x d_a
                scale = 2.0 ** (-d)
                P = G / scale                        # should be an orthogonal projection
                # projection test: P^2 = P, P = P^*, eigenvalues in {0,1}
                ev = np.linalg.eigvalsh(P)
                proj_err = np.max(np.abs(P @ P - P))
                herm_err = np.max(np.abs(P - P.conj().T))
                rk = int(np.sum(ev > 0.5))
                da, db = Va.shape[1], Vb.shape[1]
                print(f"    (a={a},b={b}) d={d}: d_a={da} d_b={db}  "
                      f"BB*/2^-d is projection? ||P^2-P||={proj_err:.1e} ||P-P*||={herm_err:.1e}  "
                      f"eig in [{ev.min():.3f},{ev.max():.3f}] rank(P)={rk} (=d_a? {rk==da})")

    print()
    print("=" * 100)
    print("(2) LOWER/DIAGONAL back-flow is ENTIRELY the r* defect:  P_a U_clean P_b = 0 for a>=b ?")
    print("=" * 100)
    for k in [8, 10, 12]:
        Uc = U_clean_matrix(k)
        bs, mz = bases_dict(k)
        n = len(mz)
        maxlow = 0.0
        for ja, a in enumerate(mz):
            for jb, b in enumerate(mz):
                if a >= b:
                    B = bs[a].conj().T @ (Uc @ bs[b])
                    maxlow = max(maxlow, np.linalg.norm(B, 2))
        # full op lower-triangular norm
        Uf = build_T(k).T
        Qfull = np.zeros((n, n))
        for ja, a in enumerate(mz):
            for jb, b in enumerate(mz):
                Qfull[ja, jb] = np.linalg.norm(bs[a].conj().T @ (Uf @ bs[b]), 2)
        lowfull = np.linalg.norm(np.tril(Qfull), 2)
        print(f"  k={k}:  max ||P_a U_CLEAN P_b||_2 over a>=b = {maxlow:.2e}  (==0 => clean op is strictly upper)"
              f"   |  ||tril(Q_full)||_2 = {lowfull:.6f} = {lowfull*2**(k/2):.4f} 2^-k/2")

    print()
    print("=" * 100)
    print("(3) ENERGY / NORM reconciliation -- closed arithmetic for the CLEAN upper block")
    print("    Q_mass[a,b] = (1/d_b)||B||_F^2,  B rank=d_b, all sv = 2^{-d/2}")
    print("    => ||B||_F^2 = d_b * 2^{-d},  Q_mass = 2^{-d}  (EXACT),  ||B||_2 = 2^{-d/2}  (EXACT).")
    print("    Dimensions: level m has d_m = 2^{k-2-m}.  rank(B_{a<b}) = d_b (full source).")
    print("=" * 100)
    for k in [8, 10, 12]:
        Uc = U_clean_matrix(k)
        bs, mz = bases_dict(k)
        print(f"\n  k={k}:  (a,b) d   d_a=2^(k-2-a) d_b=2^(k-2-b)  rank(B)  ||B||_F^2   d_b*2^-d   Q_mass   2^-d   ||B||_2  2^-d/2")
        for (a, b) in [(0, 1), (0, 3), (1, 4), (2, 5)]:
            if a in bs and b in bs and a < b:
                Va, Vb = bs[a], bs[b]
                B = Va.conj().T @ (Uc @ Vb)
                d = b - a
                da, db = Va.shape[1], Vb.shape[1]
                sv = np.linalg.svd(B, compute_uv=False)
                rk = int(np.sum(sv > 1e-9 * sv[0]))
                fro2 = np.linalg.norm(B, "fro") ** 2
                qmass = fro2 / db
                print(f"        ({a},{b}) {d}   {da:6d}       {db:6d}        {rk:5d}    {fro2:9.4f}  "
                      f"{db*2.0**-d:9.4f}  {qmass:.5f}  {2.0**-d:.5f}  {sv[0]:.5f}  {2.0**(-d/2):.5f}")

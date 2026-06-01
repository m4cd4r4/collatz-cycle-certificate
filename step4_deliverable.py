# -*- coding: utf-8 -*-
"""
STEP-4 DELIVERABLE TABLES (consolidated, all validated against build_T).
Produces exactly the three numerical artifacts the task asks for, plus the certificate row-sum.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from analytic_proofs import build_T, v2

np.set_printoptions(suppress=True, precision=5, linewidth=200)


def odds(k):
    return np.array([2 * s + 1 for s in range(1 << (k - 1))], dtype=np.int64)


def char_matrix(k):
    mod = 1 << k; N = mod >> 1
    w = np.exp(2j * np.pi / mod)
    return (w ** (np.outer(odds(k), np.arange(N)))) / np.sqrt(N)


def level_of_xi(k):
    N = 1 << (k - 1)
    return np.array([(v2(int(x)) if x != 0 else -1) for x in range(N)])


def Q_full(k):
    U = build_T(k).T
    C = char_matrix(k); lev = level_of_xi(k)
    mz = sorted(set(int(l) for l in lev if l >= 0))
    bs = {m: C[:, lev == m] for m in mz}; n = len(mz)
    Qn = np.zeros((n, n)); R = np.zeros((n, n), dtype=int); Qm = np.zeros((n, n))
    for jb, b in enumerate(mz):
        db = bs[b].shape[1]
        UVb = U @ bs[b]
        for ja, a in enumerate(mz):
            B = bs[a].conj().T @ UVb
            sv = np.linalg.svd(B, compute_uv=False)
            Qn[ja, jb] = sv[0] if sv.size else 0
            Qm[ja, jb] = np.sum(sv ** 2) / db
            R[ja, jb] = int(np.sum(sv > 1e-9 * max(sv[0], 1e-300)))
    return Qn, Qm, R, mz


def rowsum_bound(Q, t=2.0):
    n = Q.shape[0]; D = t ** np.arange(n)
    return float(np.max(((D[:, None] * Q) / D[None, :]).sum(axis=1)))


print("=" * 95)
print("TABLE 1.  Q[a,b]=||P_a U P_b||_2  vs  2^{-(b-a)/2} (UPPER a<b)  and  2^{-k/2} (LOWER+DIAG a>=b)")
print("=" * 95)
print(f"{'k':>3} | {'max|Q_up - 2^(-d/2)|':>20} | {'||tril Q||_2':>13} | {'/2^(-k/2)':>10} | {'rho(Q)':>7} | {'rowsum(t=2)':>11}")
for k in range(6, 13):
    Qn, Qm, R, mz = Q_full(k)
    n = len(mz)
    up = max((abs(Qn[a, b] - 2 ** (-(b - a) / 2)) for a in range(n) for b in range(a + 1, n)), default=0)
    low = np.linalg.norm(np.tril(Qn), 2)
    rho = float(np.max(np.abs(np.linalg.eigvals(Qn))))
    rs = rowsum_bound(Qn, 2.0)
    print(f"{k:>3} | {up:>20.3e} | {low:>13.6e} | {low*2**(k/2):>10.4f} | {rho:>7.4f} | {rs:>11.4f}")

print()
print("=" * 95)
print("TABLE 2.  RANK PROFILE of P_a U P_b.  UPPER (a<b): rank=d_b (FULL source) => NOT rank-1.")
print("          The block is 2^{-d/2} * (isometry level b -> level a).  k=10 shown.")
print("=" * 95)
Qn, Qm, R, mz = Q_full(10)
n = len(mz)
print("  rank(P_a U P_b), rows a=target level, cols b=source level:")
print(R)
print(f"\n  d_b (source dim per level) = {[2**(10-2-b) for b in mz]}")
print("  upper blocks (a<b): rank == d_b ?", all(R[a, b] == 2**(10-2-mz[b]) for a in range(n) for b in range(a+1, n)))

print()
print("=" * 95)
print("TABLE 3.  ENERGY <-> NORM reconciliation.  Q_mass[a,b]=2^{-d} EXACT, Q_norm=2^{-d/2} EXACT.")
print("          For the rank=d_b isometry block:  ||B||_F^2 = d_b*2^{-d}, Q_mass=2^{-d}, ||B||_2=2^{-d/2}.")
print("=" * 95)
print(f"  {'(a,b)':>7} {'d':>2} {'Q_mass':>9} {'2^-d':>9} {'Q_norm':>9} {'2^-d/2':>9} {'sqrt(Qmass)':>11}")
Qn, Qm, R, mz = Q_full(12)
for (a, b) in [(0, 1), (0, 2), (0, 3), (1, 3), (2, 5)]:
    ia, ib = mz.index(a), mz.index(b); d = b - a
    print(f"  {str((a,b)):>7} {d:>2} {Qm[ia,ib]:>9.5f} {2**-d:>9.5f} {Qn[ia,ib]:>9.5f} {2**(-d/2):>9.5f} {np.sqrt(Qm[ia,ib]):>11.5f}")
print("\n  Reconciliation: Q_mass is the PER-INPUT energy (Frobenius^2 / d_b); the energy spreads over")
print("  rank=d_b orthogonal output directions each carrying amplitude 2^{-d/2}.  The OPERATOR norm")
print("  picks the single worst direction = 2^{-d/2}, NOT sqrt of the summed energy.  No contradiction:")
print("  Q_mass = (avg over inputs of output energy) = 2^{-d};  Q_norm = (max amplitude) = 2^{-d/2}.")

print()
print("=" * 95)
print("CERTIFICATE row-sum  max_a sum_b Q[a,b] 2^{a-b}  (Step-3 bound on rho(Q)); must be < 1, flat in k")
print("=" * 95)
for k in range(6, 14):
    Qn, Qm, R, mz = Q_full(k)
    print(f"  k={k:2d}: row-sum(D=diag 2^j) = {rowsum_bound(Qn, 2.0):.4f}   rho(Q) = {float(np.max(np.abs(np.linalg.eigvals(Qn)))):.4f}")

# -*- coding: utf-8 -*-
"""
Pin the EXACT mechanism behind the two lemmas, so the proof route is constructive (not just observed).

LEMMA A mechanism (upper block = scaled partial isometry):
  Recall g_eta(r) = [v(r)<k][v(r)<=b] w^{eta q(r)},  q(r) = (3r+1)/2^{v(r)},  b=v2(eta).
  The map r -> q(r) is the "deterministic x3 + carry" map.  Claim to verify:
    For two source characters eta, eta' of the SAME level b, the level-a Fourier coefficients of
    g_eta and g_eta' are related by a UNITARY reindex (the x3-unit permutation), and the level-a
    energy of g_eta is EXACTLY 2^{-d} per unit input.  Concretely we verify the Gram identity
       <P_a g_eta, P_a g_eta'> = 2^{-d} delta(reindex(eta), eta')      (d=b-a)
    which is precisely  B^* B = 2^{-d} * (permutation-projection) => partial isometry of scale 2^{-d/2}.

LEMMA B mechanism (r* defect column):
  The defect operator D = U_full - U_clean is supported on the single row r=r*.  D = e_{r*} c^*,
  rank 1, where c is the "defect covector": c_eta = (U chi_eta)(r*) - g_eta(r*).
  ||D||_2 = ||c||_2 = (the L2 norm of that single row).  Verify ||c||_2 = beta 2^{-k/2}, and that the
  level-matrix of D is rank 1 with the same norm (so tril(Q_full) <= ||D in level basis|| = beta 2^{-k/2}).
  Identify beta and its k-limit.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from analytic_proofs import build_T, v2, syr

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


if __name__ == "__main__":
    print("=" * 100)
    print("LEMMA A mechanism:  B^* B = 2^{-d} * (orthogonal projection), the partial-isometry identity")
    print("  verifying it via the within-level Gram matrix of {P_a g_eta : eta in level b}")
    print("=" * 100)
    for k in [8, 10, 12]:
        Uc = U_clean_matrix(k)
        C = char_matrix(k); lev = level_of_xi(k)
        mz = sorted(set(int(l) for l in lev if l >= 0))
        bs = {m: C[:, lev == m] for m in mz}
        print(f"\n  k={k}:")
        for (a, b) in [(0, 1), (0, 2), (1, 4)]:
            if a in bs and b in bs and a < b:
                Va, Vb = bs[a], bs[b]
                B = Va.conj().T @ (Uc @ Vb)        # d_a x d_b
                d = b - a
                BstarB = B.conj().T @ B            # d_b x d_b
                scale = 2.0 ** (-d)
                Pi = BstarB / scale
                ev = np.linalg.eigvalsh(Pi)
                proj_err = np.max(np.abs(Pi @ Pi - Pi))
                db = Vb.shape[1]
                rk = int(np.sum(ev > 0.5))
                # is Pi a PERMUTATION-like 0/1 structure (each source maps to one target)?
                offdiag = Pi - np.diag(np.diag(Pi))
                print(f"    (a={a},b={b}) d={d}:  B*B/2^-d projection? ||Pi^2-Pi||={proj_err:.1e}  "
                      f"rank(Pi)={rk} of d_b={db}  diag(Pi) in [{np.diag(Pi).real.min():.3f},{np.diag(Pi).real.max():.3f}]  "
                      f"max|offdiag|={np.max(np.abs(offdiag)):.1e}")
                # so for a<b the projection Pi has rank d_b = FULL => B*B = 2^-d I_{d_b} => B is exactly
                # 2^{-d/2} times an isometry of level b into level a.
                isom_err = np.max(np.abs(BstarB - scale * np.eye(db)))
                print(f"               => B*B = 2^-d * I_{{d_b}} ?  max|B*B - 2^-d I| = {isom_err:.2e}  "
                      f"(if ~0: B is EXACTLY 2^-d/2 times an isometry; ||B||_2 = 2^-d/2 EXACT)")

    print()
    print("=" * 100)
    print("LEMMA B mechanism:  defect D = U_full - U_clean = e_{r*} c^*  (rank 1, single row).")
    print("  ||D||_2 = ||c||_2 = ||row r* of (U_full - U_clean)||_2 = beta * 2^{-k/2}.")
    print("=" * 100)
    for k in range(6, 14):
        Uf = build_T(k).T
        Uc = U_clean_matrix(k)
        D = Uf - Uc
        o = odds(k)
        v = np.array([v2(3 * int(r) + 1) for r in o])
        rstar_idx = int(np.where(v >= k)[0][0])
        rstar = int(o[rstar_idx])
        crow = D[rstar_idx, :]                  # the single nonzero row
        c_l2 = np.linalg.norm(crow)
        # level-basis level-matrix of D and its 2-norm
        C = char_matrix(k); lev = level_of_xi(k)
        mz = sorted(set(int(l) for l in lev if l >= 0))
        bs = {m: C[:, lev == m] for m in mz}; n = len(mz)
        QD = np.zeros((n, n))
        for ja, a in enumerate(mz):
            for jb, bb in enumerate(mz):
                QD[ja, jb] = np.linalg.norm(bs[a].conj().T @ (D @ bs[bb]), 2)
        svQD = np.linalg.svd(QD, compute_uv=False)
        print(f"  k={k:2d}: r*={rstar:5d}  ||c||_2={c_l2:.6f}={c_l2*2**(k/2):.4f}*2^-k/2  "
              f"rank(D)={np.linalg.matrix_rank(D,tol=1e-9)}  "
              f"||Q_D||_2={svQD[0]:.6f}={svQD[0]*2**(k/2):.4f}*2^-k/2 rank(Q_D)={int(np.sum(svQD>1e-9*svQD[0]))}")

    print()
    print("=" * 100)
    print("Lemma B closed form attempt:  the r* row is (U chi_eta)(r*) = (1/2^k) sum_m w^{eta Syr(r*+m2^k)}")
    print("  Its L2-over-eta norm.  Check ||c||_2^2 = (1/2^k) sum_{nonzero-level eta} |..|^2 -> beta^2 2^-k.")
    print("=" * 100)
    for k in [8, 10, 12]:
        Uf = build_T(k).T
        Uc = U_clean_matrix(k)
        D = Uf - Uc
        o = odds(k); v = np.array([v2(3*int(r)+1) for r in o])
        ri = int(np.where(v >= k)[0][0])
        crow = D[ri, :]
        # split contribution by source level b = v2(eta)
        lev_eta = np.array([v2(e) if e != 0 else -1 for e in range(len(crow))])
        print(f"\n  k={k}: defect row energy by source level b:")
        tot = 0.0
        for b in range(0, k - 1):
            mask = lev_eta == b
            e = np.sum(np.abs(crow[mask]) ** 2)
            tot += e
            if e > 1e-12:
                print(f"      b={b}: energy={e:.6e}  (x2^k = {e*2**k:.4f})")
        print(f"      total ||c||^2={tot:.6e}  ||c||={np.sqrt(tot):.6f}  beta={np.sqrt(tot)*2**(k/2):.4f}")

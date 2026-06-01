# -*- coding: utf-8 -*-
"""
STEP-4 FOUNDATION (corrected).  Exact block-formula derivation of  Q[a,b] = ||P_a U P_b||_2.

KEY CORRECTION to the prompt's stated formula
   (U chi_eta)(r) = [v(r)<=b] w^{eta(3r+1)/2^{v(r)}}        <-- WRONG as literally stated
The honest exact identity (verified to 3e-16 against build_T) is:
   (U chi_eta)(r) = [v(r)<=b] w^{eta(3r+1)/2^{v(r)}}   ONLY for v(r) < k ;
   for v(r) >= k (the single exceptional residue r*),  (U chi_eta)(r*) is a genuine partial
   Gauss sum, NOT given by the masked-phase formula.  That single defect IS the lower-block back-flow.

This file proves/derives, and validates in Python:
  (A) the corrected exact block formula (already cross-checked in diagnose_block_formula.py);
  (B) the level decomposition of  g_eta(r) := [v(r)<k][v(r)<=b] w^{eta(3r+1)/2^{v(r)}}  and the
      level-a amplitude  ->  ||P_a (clean part) P_b||_2 has the CLOSED FORM 2^{-(b-a)/2};
  (C) the singular-value structure: upper blocks are NOT rank-1; they are SCALED PARTIAL ISOMETRIES
      (all nonzero singular values equal to 2^{-(b-a)/2}), so ||.||_2 = 2^{-(b-a)/2} exactly;
  (D) the r* defect: the lower/diagonal back-flow is exactly rank-1 of scale 2^{-k/2}.

Conventions identical to level_transfer_Q.py.  U=T_k^T,  Q[a,b]=||P_a U P_b||_2, a=target,b=source.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from analytic_proofs import build_T, v2, syr

np.set_printoptions(suppress=True, precision=6, linewidth=220)


# ------------------------- frequency machinery -------------------------
def odds(k):
    return np.array([2 * s + 1 for s in range(1 << (k - 1))], dtype=np.int64)


def char_matrix(k):
    mod = 1 << k; N = mod >> 1
    w = np.exp(2j * np.pi / mod)
    return (w ** (np.outer(odds(k), np.arange(N)))) / np.sqrt(N)


def level_of_xi(k):
    N = 1 << (k - 1)
    return np.array([(v2(int(x)) if x != 0 else -1) for x in range(N)])


# ------------------------- the CLEAN operator (drop the r* defect) -------------------------
def U_clean_on_chi(k, eta):
    """
    g_eta(r) = [v(r) < k] [v(r) <= b] w^{eta(3r+1)/2^{v(r)}},   b=v2(eta).
    This is U chi_eta with the r* row zeroed (the 'clean' coset-collapse part).
    """
    mod = 1 << k; N = mod >> 1
    w = np.exp(2j * np.pi / mod)
    b = v2(eta)
    o = odds(k)
    out = np.zeros(N, dtype=complex)
    for i in range(N):
        r = int(o[i]); v = v2(3 * r + 1)
        if v < k and v <= b:
            q = (3 * r + 1) >> v
            out[i] = w ** ((eta * q) % mod)
    return out


def U_clean_matrix(k):
    """
    Clean operator in the STATE basis = U_full with the exceptional r* rows zeroed.
    Rationale: (U f)(r) = average of f over the Syracuse coset of r.  For v(r) < k the coset is a
    FULL coset (uniform, multiplicity 2^{k-v}); the Gauss-sum collapse applies and gives the masked
    phase.  Only r with v(r) >= k (the unique r* and -- at the working modulus -- its sibling) fail
    coset-uniformity.  So  U_clean := U_full restricted to rows r with v(r) < k.
    """
    mod = 1 << k; N = mod >> 1
    o = odds(k)
    v = np.array([v2(3 * int(r) + 1) for r in o])
    Uc = Ufull_matrix(k).copy()
    Uc[v >= k, :] = 0.0          # zero the rows where coset-uniformity fails (the r* defect rows)
    return Uc


def Ufull_matrix(k):
    return build_T(k).T


# ------------------------- block / Q machinery -------------------------
def Q_blocks(k, U):
    """Given the operator U as an N x N matrix in the STATE basis, return Qnorm, Qmass, ranks,
    and per-block singular spectra in the level basis."""
    C = char_matrix(k); lev = level_of_xi(k)
    mz = sorted(set(int(l) for l in lev if l >= 0))
    bases = {m: C[:, lev == m] for m in mz}
    L = len(mz)
    Qn = np.zeros((L, L)); Qm = np.zeros((L, L)); R = np.zeros((L, L), dtype=int)
    svspec = {}
    UV = {m: U @ bases[m] for m in mz}
    for jb, b in enumerate(mz):
        db = bases[b].shape[1]
        for ja, a in enumerate(mz):
            B = bases[a].conj().T @ UV[b]
            sv = np.linalg.svd(B, compute_uv=False)
            Qn[ja, jb] = sv[0] if sv.size else 0.0
            Qm[ja, jb] = (np.linalg.norm(B, "fro") ** 2) / db
            tol = 1e-9 * max(sv[0], 1e-300) if sv.size else 1
            R[ja, jb] = int(np.sum(sv > tol))
            svspec[(a, b)] = sv
    return {"levels": mz, "Qnorm": Qn, "Qmass": Qm, "ranks": R, "svspec": svspec,
            "bases": bases}


# ============================ MAIN ============================
if __name__ == "__main__":
    print("#" * 100)
    print("# (A0) consistency: state-basis U_clean (rows v(r)>=k zeroed) == per-column g_eta formula?")
    print("#" * 100)
    for k in range(4, 11):
        N = 1 << (k - 1)
        Uc = U_clean_matrix(k)
        C = char_matrix(k)
        max_err = 0.0
        for eta in range(1, N):
            chi = C[:, eta] * np.sqrt(N)          # un-normalized chi_eta(r)=w^{eta r}
            lhs = Uc @ chi
            rhs = U_clean_on_chi(k, eta)
            max_err = max(max_err, np.max(np.abs(lhs - rhs)))
        o = odds(k)
        v = np.array([v2(3 * int(r) + 1) for r in o])
        rstar = [int(o[i]) for i in np.where(v >= k)[0]]
        print(f"  k={k:2d}: max|U_clean@chi - g_eta| = {max_err:.2e}   #rows zeroed (v>=k) = {int(np.sum(v>=k))}  r*={rstar}")

    print()
    print("#" * 100)
    print("# (A) CORRECTED exact block formula:  U = U_clean  +  defect(r*).  Verify decomposition.")
    print("#" * 100)
    for k in range(4, 12):
        Ufull = Ufull_matrix(k)
        Uclean = U_clean_matrix(k)
        D = Ufull - Uclean                    # the defect operator
        # defect should be supported on the r* rows only
        rownorm = np.linalg.norm(D, axis=1)
        nz_rows = np.where(rownorm > 1e-12)[0]
        o = odds(k)
        rstar_rows = [int(o[i]) for i in nz_rows]
        # rank of defect
        dr = np.linalg.matrix_rank(D, tol=1e-9)
        print(f"  k={k:2d}: ||U_full - U_clean||_2 = {np.linalg.norm(D,2):.6e} = "
              f"{np.linalg.norm(D,2)*2**(k/2):.4f}*2^(-k/2) | #defect nz-rows={len(rstar_rows)} r={rstar_rows} | rank(defect)={dr}")

    print()
    print("#" * 100)
    print("# (B) CLEAN operator U_clean: each block P_a U_clean P_b  -- closed form 2^{-(b-a)/2}?")
    print("#     and the singular-value structure (rank-1 ?  or scaled partial isometry ?)")
    print("#" * 100)
    for k in [8, 10, 12]:
        bl = Q_blocks(k, U_clean_matrix(k))
        Qn = bl["Qnorm"]; R = bl["ranks"]; mz = bl["levels"]; n = len(mz)
        # upper deviation from 2^{-d/2}
        dev = max((abs(Qn[a, b] - 2 ** (-(b - a) / 2)) for a in range(n) for b in range(a + 1, n)),
                  default=0.0)
        # singular-value flatness of a representative upper block: are ALL nonzero sv equal?
        a, b = 0, 2
        sv = bl["svspec"][(a, b)]
        sv_nz = sv[sv > 1e-9 * sv[0]]
        flat = (sv_nz.max() - sv_nz.min())
        print(f"  k={k:2d}: clean upper-cascade max|Q[a,b]-2^(-d/2)|={dev:.2e}   "
              f"block(0,2): rank={R[0,2]} of d_b={bl['bases'][2].shape[1]}, "
              f"#nonzero-sv={sv_nz.size}, all-equal? spread={flat:.2e}, common sv={sv_nz.mean():.6f} (=2^-1={0.5})")

    print()
    print("#" * 100)
    print("# (C) Q[a,b] = ||P_a U_full P_b||_2  vs targets, k=6..12  (THE deliverable table)")
    print("#" * 100)
    print("  k |  upper a<b: max|Q-2^(-d/2)|  | lower+diag a>=b: ||tril Q||_2 (x 2^(k/2))")
    for k in range(6, 13):
        bl = Q_blocks(k, Ufull_matrix(k))
        Qn = bl["Qnorm"]; n = Qn.shape[0]
        dev = max((abs(Qn[a, b] - 2 ** (-(b - a) / 2)) for a in range(n) for b in range(a + 1, n)),
                  default=0.0)
        low = np.tril(Qn); ln = np.linalg.norm(low, 2)
        print(f"  {k:2d} |   {dev:.3e}              |  {ln:.6e}  ({ln*2**(k/2):.4f} * 2^(-k/2))")

    print()
    print("#" * 100)
    print("# (D) r*-defect / lower-diagonal back-flow block: rank and 2^{-k/2} scale")
    print("#" * 100)
    for k in [8, 10, 12]:
        bl = Q_blocks(k, Ufull_matrix(k))
        Qn = bl["Qnorm"]
        low = np.tril(Qn)
        sv = np.linalg.svd(low, compute_uv=False)
        print(f"  k={k:2d}: tril(Q) singular values [:4]={np.round(sv[:4],6)}  "
              f"||tril||_2={sv[0]:.6f}={sv[0]*2**(k/2):.4f}*2^(-k/2)  "
              f"sv[1]/sv[0]={sv[1]/sv[0]:.4f}")
        # The DEFECT itself (U_full - U_clean) in level basis
        D = Ufull_matrix(k) - U_clean_matrix(k)
        blD = Q_blocks(k, D)
        QnD = blD["Qnorm"]
        svD = np.linalg.svd(QnD, compute_uv=False)
        print(f"        defect-operator level matrix Qnorm: sv[:3]={np.round(svD[:3],6)}  "
              f"rank={int(np.sum(svD>1e-9*svD[0]))}  ||.||_2*2^(k/2)={svD[0]*2**(k/2):.4f}")

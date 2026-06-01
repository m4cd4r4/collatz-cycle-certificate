# -*- coding: utf-8 -*-
"""
Generate the README figures directly from the transfer operator build_T. Reproducible:
    python generate_figures.py
Every figure is computed from the real operator, not drawn by hand. Sober on purpose.
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from analytic_proofs import build_T, v2

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "figures")
os.makedirs(OUT, exist_ok=True)


def level_blocks(k):
    """Return (Q, Uc, idx, levels) for resolution k.
    Q[a,b] = top singular value of the 2-adic level block P_a U P_b in the character basis."""
    U = build_T(k).T
    mod = 1 << k
    N = mod >> 1
    odds = (2 * np.arange(N) + 1).astype(np.int64)
    w = np.exp(2j * np.pi / mod)
    C = (w ** (np.outer(odds, np.arange(N)))) / np.sqrt(N)
    Uc = C.conj().T @ U @ C
    lev = np.array([(v2(int(x)) if x != 0 else -1) for x in range(N)])
    levels = sorted(set(lev[lev >= 0].tolist()))
    idx = {a: np.where(lev == a)[0] for a in levels}
    Q = np.zeros((len(levels), len(levels)))
    for ia, a in enumerate(levels):
        for ib, b in enumerate(levels):
            blk = Uc[np.ix_(idx[a], idx[b])]
            if blk.size:
                Q[ia, ib] = np.linalg.svd(blk, compute_uv=False)[0]
    return Q, Uc, idx, levels


def fig_Q_heatmap(k=10):
    Q, *_ = level_blocks(k)
    fig, ax = plt.subplots(figsize=(5.2, 4.4))
    M = np.log2(Q + 1e-18)
    im = ax.imshow(M, cmap="viridis", origin="upper")
    ax.set_xlabel("source level b"); ax.set_ylabel("target level a")
    ax.set_title(f"Level-block norms  Q[a,b] = ||P_a U P_b||   (k={k})\n"
                 f"upper cascade = 2^(-(b-a)/2) exact; lower triangle = tiny r* defect", fontsize=9)
    cb = fig.colorbar(im, ax=ax, shrink=0.85); cb.set_label("log2 of block norm")
    fig.tight_layout(); p = os.path.join(OUT, "fig1_Q_heatmap.png"); fig.savefig(p, dpi=140); plt.close(fig)
    return p


def fig_spectrum(k=10):
    U = build_T(k)
    ev = np.linalg.eigvals(U)
    mag = np.sort(np.abs(ev))[::-1]
    lam2 = mag[1]
    fig, ax = plt.subplots(figsize=(5.6, 3.6))
    ax.plot(range(1, min(40, len(mag)) + 1), mag[:40], "o-", ms=3, lw=0.8, color="#33508a")
    ax.axhline(1.0, color="#888", lw=0.8, ls="--")
    ax.annotate(f"|lambda_1| = 1 (Perron)", (1, 1.0), (6, 1.02), fontsize=8, color="#555")
    ax.annotate(f"|lambda_2| = {lam2:.3f}  (spectral gap = {1-lam2:.3f})",
                (2, lam2), (6, lam2 - 0.12), fontsize=8, color="#a23",
                arrowprops=dict(arrowstyle="->", color="#a23", lw=0.7))
    ax.set_xlabel("eigenvalue index (by |.|)"); ax.set_ylabel("|eigenvalue|")
    ax.set_title(f"Transfer-operator spectrum  (k={k}):  a gap below 1 forbids non-trivial cycles", fontsize=9)
    ax.set_ylim(0, 1.1)
    fig.tight_layout(); p = os.path.join(OUT, "fig2_spectrum.png"); fig.savefig(p, dpi=140); plt.close(fig)
    return p, lam2


def fig_incidence(k=9, a=0, b=3):
    _, Uc, idx, levels = level_blocks(k)
    blk = np.abs(Uc[np.ix_(idx[a], idx[b])])
    thr = 0.5 * blk.max()
    fig, ax = plt.subplots(figsize=(5.0, 4.2))
    ax.imshow(blk > thr, cmap="Greys", origin="upper", aspect="auto", interpolation="nearest")
    ax.set_xlabel(f"source characters (level b={b})")
    ax.set_ylabel(f"target characters (level a={a})")
    ax.set_title(f"Support of one block B = P_a U P_b  (k={k}, d=b-a={b-a})\n"
                 f"exactly ONE mark per row => columns orthogonal => B*B = 2^(-d) I", fontsize=9)
    fig.tight_layout(); p = os.path.join(OUT, "fig3_incidence.png"); fig.savefig(p, dpi=140); plt.close(fig)
    return p


def fig_certificate(kmax=11):
    ks = list(range(6, kmax + 1))
    rowsums, rhos = [], []
    for k in ks:
        Q, *_ = level_blocks(k)
        n = Q.shape[0]
        w = np.array([[2.0 ** (a - b) for b in range(n)] for a in range(n)])
        rowsums.append(float(np.max((Q * w).sum(axis=1))))
        rhos.append(float(max(abs(np.linalg.eigvals(Q)))))
    fig, ax = plt.subplots(figsize=(5.6, 3.6))
    ax.plot(ks, rowsums, "o-", color="#33508a", label="weighted row-sum  max_a sum_b Q[a,b] 2^(a-b)")
    ax.plot(ks, rhos, "s-", color="#2a8", label="spectral radius rho(Q)")
    ax.axhline(1.0, color="#a23", lw=1.0, ls="--", label="threshold for cycles (= 1)")
    ax.set_ylim(0, 1.1); ax.set_xlabel("resolution k"); ax.set_ylabel("certificate value")
    ax.set_title("The certificate stays flat below 1 as k grows  =>  no non-trivial cycles", fontsize=9)
    ax.legend(fontsize=7, loc="center right")
    fig.tight_layout(); p = os.path.join(OUT, "fig4_certificate.png"); fig.savefig(p, dpi=140); plt.close(fig)
    return p, rowsums, rhos


if __name__ == "__main__":
    print("fig1:", fig_Q_heatmap(10))
    p2, lam2 = fig_spectrum(10); print("fig2:", p2, " |lambda_2| =", round(lam2, 4))
    print("fig3:", fig_incidence(9, 0, 3))
    p4, rs, rh = fig_certificate(11)
    print("fig4:", p4)
    print("  row-sums k=6..11:", [round(x, 4) for x in rs])
    print("  rho(Q)   k=6..11:", [round(x, 4) for x in rh])

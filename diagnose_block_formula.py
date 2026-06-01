# -*- coding: utf-8 -*-
"""
Diagnose the EXACT action of U = build_T(k).T on a character chi_eta, residue by residue,
to pin the correct closed form. The prompt's claimed formula
   (U chi_eta)(r) = [v(r)<=b] w^{eta(3r+1)/2^{v(r)}}
FAILS against build_T (err ~0.25). Find what build_T actually computes.

build_T column-stochastic def (analytic_proofs.build_T):
   for src in range(N): r=2*src+1
     for m in range(mod): n = r + m*mod;  s = syr(n) % mod;  tgt=(s-1)//2;  T[tgt,src]+=1
   T /= mod
So T[tgt,src] = (1/2^k) #{ m in [0,2^k) : syr(2*src+1 + m*2^k) == 2*tgt+1 (mod 2^k) }.
U = T^T:  (U f)(src) = sum_tgt T[tgt,src] f(tgt)
        = (1/2^k) sum_{m=0}^{2^k-1} f_state( (syr(r+m*2^k)-1)/2 )
i.e. (U f)(r) = average over the 2^k lifts m of f evaluated at the Syracuse image.
For f = chi_eta as a STATE function: chi_eta(state=tgt) = w^{eta*(2*tgt+1)} = w^{eta * s}, s=2tgt+1 = syr(...) mod 2^k.
So  (U chi_eta)(r) = (1/2^k) sum_{m} w^{ eta * (syr(r+m 2^k) mod 2^k) }.
THAT is the exact object. Let's compute it and compare to candidate closed forms.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from analytic_proofs import build_T, v2, syr

np.set_printoptions(suppress=True, precision=5, linewidth=200)


def exact_U_chi(k, eta):
    """Direct from definition: (U chi_eta)(r) = (1/2^k) sum_m w^{eta * syr(r+m 2^k) mod 2^k}."""
    mod = 1 << k
    N = mod >> 1
    w = np.exp(2j * np.pi / mod)
    out = np.zeros(N, dtype=complex)
    for src in range(N):
        r = 2 * src + 1
        acc = 0j
        for m in range(mod):
            s = syr(r + m * mod) % mod
            acc += w ** ((eta * s) % mod)
        out[src] = acc / mod
    return out


def via_buildT(k, eta):
    mod = 1 << k
    N = mod >> 1
    w = np.exp(2j * np.pi / mod)
    U = build_T(k).T
    chi = np.array([w ** ((eta * (2 * src + 1)) % mod) for src in range(N)], dtype=complex)
    return U @ chi


def candidate_gausssum(k, eta):
    """
    Honest Gauss-sum collapse. For r with v=v2(3r+1):
      Syr(r + m 2^k) ranges over the coset { (3r+1)/2^v + j 2^{k-v} : j=0..2^v-1 } UNIFORMLY,
      each value hit 2^{k-v} times (the coset-uniformity lemma), PROVIDED v < k.
    Then (1/2^k) sum_m w^{eta*Syr} = (1/2^v) sum_{j=0}^{2^v-1} w^{ eta ((3r+1)/2^v + j 2^{k-v}) }
       = w^{eta (3r+1)/2^v} * (1/2^v) sum_j w^{eta j 2^{k-v}}.
    The inner sum is a geometric/Gauss sum:
       (1/2^v) sum_{j=0}^{2^v-1} W^{j},  W = w^{eta 2^{k-v}} = exp(2 pi i eta 2^{k-v} / 2^k)
                                                            = exp(2 pi i eta / 2^v).
    This is (1/2^v) * [ 2^v if 2^v | eta else (W^{2^v}-1)/(W-1) ].  W^{2^v}=exp(2 pi i eta)=1 always,
    so the sum is 2^v * [2^v | eta] = 2^v * [ v2(eta) >= v ] = 2^v * [ v <= b ].
    => (U chi_eta)(r) = [ v <= b ] * w^{ eta (3r+1)/2^v }    when v<k.
    For v>=k (only r=r*): Syr is NOT spread over a full coset; handle separately (return the direct value).
    """
    mod = 1 << k
    N = mod >> 1
    w = np.exp(2j * np.pi / mod)
    b = v2(eta)
    out = np.zeros(N, dtype=complex)
    for src in range(N):
        r = 2 * src + 1
        v = v2(3 * r + 1)
        if v < k:
            if v <= b:
                q = (3 * r + 1) >> v
                out[src] = w ** ((eta * q) % mod)
            else:
                out[src] = 0.0
        else:
            # exceptional r* : compute directly
            acc = 0j
            for m in range(mod):
                s = syr(r + m * mod) % mod
                acc += w ** ((eta * s) % mod)
            out[src] = acc / mod
    return out


if __name__ == "__main__":
    print("Compare three computations of (U chi_eta)(r): build_T.T,  direct-from-def,  Gauss-sum candidate")
    print("(direct-from-def IS the definition, so build_T.T must match it exactly)")
    for k in range(4, 11):
        mod = 1 << k
        N = mod >> 1
        max_bt_def = 0.0
        max_def_cand = 0.0
        worst_eta = None
        for eta in range(1, N):
            d1 = exact_U_chi(k, eta)
            d2 = via_buildT(k, eta)
            d3 = candidate_gausssum(k, eta)
            e_bt = np.max(np.abs(d1 - d2))
            e_cand = np.max(np.abs(d1 - d3))
            if e_bt > max_bt_def:
                max_bt_def = e_bt
            if e_cand > max_def_cand:
                max_def_cand = e_cand
                worst_eta = eta
        print(f"  k={k:2d}:  ||buildT.T - def||={max_bt_def:.2e}   ||def - GaussSumCandidate||={max_def_cand:.2e}"
              f"   (worst eta={worst_eta}, b=v2={v2(worst_eta) if worst_eta else '-'})")

    # If the candidate matches the def, the prompt's formula was missing the v<k restriction
    # and/or the r* exception. Show the disagreement of the PROMPT's literal formula
    # (which used [v(r)<=b] with NO v<k carve-out and a different exponent handling).
    print("\nDetail at k=6, eta=1 (b=0): per-residue r, v(r), def value, candidate value")
    k = 6; eta = 1; mod = 1 << k; N = mod >> 1; w = np.exp(2j*np.pi/mod)
    d1 = exact_U_chi(k, eta); d3 = candidate_gausssum(k, eta)
    for src in range(N):
        r = 2*src+1; v = v2(3*r+1)
        if abs(d1[src]-d3[src]) > 1e-9 or v >= k:
            print(f"   r={r:3d} v={v}  def={d1[src]: .4f}  cand={d3[src]: .4f}  |diff|={abs(d1[src]-d3[src]):.2e}")

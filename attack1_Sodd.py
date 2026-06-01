# -*- coding: utf-8 -*-
"""
Closed form Sodd(alpha,m) = sum_{q odd, 0<=q<2^m} exp(2pi i alpha q/2^k),  w=exp(2pi i/2^k).
q=2l+1, l=0..2^{m-1}-1: Sodd = w^alpha sum_{l=0}^{2^{m-1}-1} w^{2 alpha l}=w^alpha * G(2alpha, 2^{m-1})
 where G(beta,L)=sum_{l=0}^{L-1} exp(2pi i beta l/2^k). G = L if 2^k|beta else (zeta^L-1)/(zeta-1), zeta=exp(2pi i beta/2^k).
For beta=2alpha, L=2^{m-1}: zeta^L=exp(2pi i 2alpha 2^{m-1}/2^k)=exp(2pi i alpha 2^m/2^k)=exp(2pi i alpha /2^{k-m}).
Let me just get the magnitude. Key question for orthogonality: when is Sodd(alpha,m) NONZERO and what's its phase.
Write alpha = 2^s * odd (s=v2(alpha)), or alpha=0.
 - The geometric sum over odd q of w^{alpha q}: standard. Let me tabulate |Sodd| and identify the rule.
"""
import numpy as np
def Sodd(alpha,m,k):
    mod=1<<k
    return sum(np.exp(2j*np.pi*(alpha%mod)*q/mod) for q in range(1<<m) if q%2==1)
def v2(n):
    if n==0: return 99
    c=0
    while n%2==0: n//=2; c+=1
    return c

k=10
print("Sodd(alpha,m): tabulate by v2(alpha) and relation of m to k.")
print("Claim to find: Sodd(alpha,m) nonzero only for specific v2(alpha); |Sodd|=2^{m-1} when alpha≡0 mod 2^k...")
for m in [4,6,8]:
    print(f" m={m}:")
    for alpha in range(0, 1<<k, 1):
        S=Sodd(alpha,m,k)
        if abs(S)>1e-9:
            s=v2(alpha)
            # record (alpha, v2, |S|, phase)
            pass
    # systematically: for each s=v2(alpha) from 0..k, pick one alpha and show |S|
    for s in range(0,k+1):
        alpha=(1<<s) if s<k else 0
        S=Sodd(alpha,m,k)
        print(f"    v2(alpha)={s:2d} (alpha={alpha:4d}): |S|={abs(S):8.4f}  S={S:.4f}")

"""Bisect the zeros of qt(1) and qt(t) in omega/m, and compare them with the CTdC endpoints."""
import mpmath as mp
from step5_window import dual_q_signs, ctdc_window

def zero_of(idx, Lam, a, M, lo, hi, s=-2, m=2, lam=7.3):
    f = lambda x: dual_q_signs(Lam, a, M, s, m*x, m, lam)[idx].real
    fa, fb = f(lo), f(hi)
    assert fa*fb < 0
    for _ in range(45):
        mid = (lo + hi)/2
        fm = f(mid)
        if fa*fm <= 0: hi, fb = mid, fm
        else: lo, fa = mid, fm
    return (lo + hi)/2

for (Lam, a, M) in [(0.05, 0.6, 1.0), (0.02, 0.9, 1.0), (0.1, 0.3, 1.0)]:
    lo, hi, Om2, Om1 = ctdc_window(Lam, a, M)
    z1 = zero_of(0, Lam, a, M, Om2*0.5, Om1*1.3)
    zt = zero_of(1, Lam, a, M, Om2*0.5, Om1*1.3)
    print(f'Lam={Lam} a={a}: zeros of qt(1), qt(t) at omega/m = {z1:.10f}, {zt:.10f};  CTdC endpoints {lo:.10f}, {hi:.10f}')

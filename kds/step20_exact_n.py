"""Step 20: n0, n1, n2 from the CLOSED-FORM leading dual P2 (step 15), with no m-extrapolation.
q(1), q(t) come from the dual horizon exponents: q(z_h)^2 = P2(z_h), with the sign pattern q(1) < 0 < q(t)."""
import warnings; warnings.filterwarnings('ignore')
import io, contextlib, numpy as np
with contextlib.redirect_stdout(io.StringIO()):
    from step15_symbolic_leading import P2_closed_form
    from step17_ellmin import ell_min
from step5_window import ctdc_window
P = np.polynomial.polynomial

def n_coeffs(Lam, a, w, ell):
    t, P2, err, info = P2_closed_form(Lam, a, 1.0, w, ell, branch=1)
    p2 = np.array([c.real for c in P2])
    q1v = -np.sqrt(max(P.polyval(1.0, p2), 0)); qtv = np.sqrt(max(P.polyval(t, p2), 0))
    s1 = (qtv - q1v)/(t - 1); q = np.array([q1v - s1, s1])
    N, rem = P.polydiv(P.polysub(P.polymul(q, q), p2), P.polymul([-1.0, 1.0], [t, -1.0]))
    N = list(N) + [0, 0, 0]
    return N[0], N[1], N[2], float(np.max(np.abs(rem))), t, P.polyval(1.0, p2), P.polyval(t, p2)

worst_n0 = worst_n1 = np.inf; maxn2 = 0
for (Lam, a) in [(0.05, 0.6), (0.02, 0.9), (0.1, 0.3), (0.001, 0.95), (0.05, 0.9877), (0.08, 0.5), (0.01, 0.99), (0.0001, 0.99), (0.001, 0.999), (0.1, 0.05), (0.11, 0.01)]:
    try:
        lo, hi, _, _ = ctdc_window(Lam, a, 1.0)
    except Exception as ex:
        print('skip', Lam, a, ex); continue
    for wf in np.linspace(0.001, 0.999, 15):
        w = lo + wf*(hi - lo); lm = ell_min(Lam, a, w)
        n0, n1, n2, rem, t, P21, P2t = n_coeffs(Lam, a, w, lm)
        scale = abs(n0) + abs(n1)
        maxn2 = max(maxn2, abs(n2)/scale)
        worst_n0 = min(worst_n0, n0/scale); worst_n1 = min(worst_n1, n1/scale)
        if P21 < 0 or P2t < 0: print('P2 negative at a horizon!', Lam, a, wf)
    print(f'Lam={Lam} a={a}: done')
print(f'max |n2|/(|n0|+|n1|) = {maxn2:.2e};  min n0/scale = {worst_n0:+.3e};  min n1/scale at ell_min = {worst_n1:+.3e}')

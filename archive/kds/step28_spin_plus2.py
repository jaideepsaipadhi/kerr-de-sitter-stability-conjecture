"""Step 28: spin s = +2.  (a) Is the (0,0,0)A dual still real?  (b) Does its superradiant window still match CTdC (3.37)?
(c) Does the explicit layer construction (step 25) close for s = +2?"""
import warnings; warnings.filterwarnings('ignore')
import io, contextlib
with contextlib.redirect_stdout(io.StringIO()):
    from step17_ellmin import ell_min
import numpy as np
from step3_dual import scan
from step5_window import dual_q_signs, ctdc_window
import step25_layer as S25

for (Lam, a) in [(0.05, 0.6), (0.02, 0.9), (0.001, 0.95)]:
    lo, hi, Om2, Om1 = ctdc_window(Lam, a, 1.0)
    res, _ = scan(Lam, a, 1.0, +2, complex(2*(lo + hi)/2), 2, 7.3)
    real = [(c, 'AB'[w]) for (c, w, imQ, reQ, imo, le, zr) in res if imQ < 1e-6]
    print(f'Lam={Lam} a={a}: real duals for s=+2: {real}')
    def zero(idx):
        f = lambda x: dual_q_signs(Lam, a, 1.0, +2, 2*x, 2, 7.3)[idx].real
        A_, B_ = Om2*0.5, Om1*1.3
        fa = f(A_)
        for _ in range(45):
            M_ = (A_ + B_)/2; fm = f(M_)
            if fa*fm <= 0: B_ = M_
            else: A_, fa = M_, fm
        return (A_ + B_)/2
    print(f'   dual window endpoints (s=+2): {sorted([zero(0), zero(1)])}  vs CTdC ({lo:.10f}, {hi:.10f})')
    oks = []
    for wf in (0.05, 0.5, 0.95):
        w = lo + wf*(hi - lo)
        for m in (1, 4, 32):
            ell = ell_min(Lam, a, w)
            ok = False
            for (d, p) in ((0.5, 0.5), (1.0, 0.3), (0.5, 0.1)):
                F = S25.build(Lam, a, m, ell, wf, s=+2)
                if S25.check(F, d, p)['ok']: ok = True; break
            oks.append(ok)
    print(f'   layer construction closes (s=+2): {sum(oks)}/{len(oks)}')

"""Certify a few physical parameter points rigorously, and confirm the verifier rejects a known-bad multiplier."""
import time
from verifier import certify_point
cases = [  # (t, z0-fraction, A, B-factor over the floor 4t - 1/2)
    (2.5, 0.5, 1.0, 1.3), (44.0, 0.2, 70.0, 1.3), (1.33, 0.9, 0.05, 1.5), (7.0, 0.95, 10.0, 1.3), (265.0, 0.5, 1e4, 1.3)]
for (t, zf, A, Bf) in cases:
    z0 = 1 + zf*(t - 1); B = Bf*(4*t - 0.5)
    L1 = A + B - 3
    pst = 2*(z0 - 1)*L1/(t - 1)
    done = False
    for c in (0.5, 0.2, 0.05):
        p = min(0.9, c*pst)
        for d in (0.5, 1.0, 0.25, 2.0):
            for b in (1.0, 0.25, 4.0, 25.0):
                t0 = time.time()
                ok, info = certify_point(t, z0, A, B, d, p, b)
                if ok:
                    print(f't={t} zf={zf} A={A} B={B:.2f}: CERTIFIED with delta={d} p={p:.4g} beta={b}  ({info} boxes, {time.time()-t0:.1f}s)')
                    done = True; break
            if done: break
        if done: break
    if not done: print(f't={t} zf={zf} A={A} B={B:.2f}: not certified')
# negative control: p > 1 violates C at z = 1, so the verifier must refuse
ok, info = certify_point(2.5, 1.75, 1.0, 12.35, 0.5, 1.2, 1.0, max_depth=18)
print('negative control (p = 1.2 must fail):', ok, info[0] if not ok else info)

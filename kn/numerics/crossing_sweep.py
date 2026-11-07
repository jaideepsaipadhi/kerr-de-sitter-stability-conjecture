# A5 decisive check: look for a coupled KN mode crossing into Im w > 0 (sign change of Im w of least-damped mode).
# Grid: rho = sqrt(a^2+Q^2) in RHOS, angle psi (a = rho cos psi, Q = rho sin psi), m in MS.
# Each point is solved at two resolutions; only eigenvalues present in BOTH (|dw| < TOL) are trusted.
# Output: crossing_sweep.csv  (a, Q, m, Re w, Im w, |dw| between resolutions) for the 3 least-damped trusted modes,
#         and a FLAG line for any trusted mode with Im w > -1e-4.
import numpy as np, sys, time, csv
from qnm_coupled import qnms
RHOS = [0.90, 0.98]
PSIS = np.linspace(0.2, np.pi/2 - 0.1, 4)      # a > 0, Q > 0
MS = [1, 2]
RES = [(20, 8), (28, 10)]
TOL = 2e-3
out = open('crossing_sweep.csv', 'w', newline=''); W = csv.writer(out)
W.writerow(['a', 'Q', 'm', 'Re_w', 'Im_w', 'res_diff'])
t0 = time.time()
for rho in RHOS:
    for psi in PSIS:
        a, Q = rho*np.cos(psi), rho*np.sin(psi)
        for m in MS:
            ev = [np.array(qnms(a, Q, m, Nr, Nt)) for (Nr, Nt) in RES]
            trusted = []
            for z in ev[1]:
                d = np.min(np.abs(ev[0] - z)) if len(ev[0]) else 1e9
                if d < TOL: trusted.append((z, d))
            trusted.sort(key=lambda t: -t[0].imag)
            for z, d in trusted[:3]:
                W.writerow([f'{a:.4f}', f'{Q:.4f}', m, f'{z.real:.6f}', f'{z.imag:.6f}', f'{d:.1e}'])
                if z.imag > -1e-4:
                    print(f'FLAG a={a:.4f} Q={Q:.4f} m={m} w={z:.6f} diff={d:.1e}', flush=True)
            out.flush()
            lead = trusted[0][0] if trusted else None
            print(f'[{time.time()-t0:7.0f}s] a={a:.3f} Q={Q:.3f} m={m} least-damped={lead}', flush=True)
print('done')

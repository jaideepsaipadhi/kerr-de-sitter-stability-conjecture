"""Margins of the large-m conditions, including toward the window edges and as ell -> infinity."""
import warnings; warnings.filterwarnings('ignore')
from step13_largem import check
from step9_freeq import valid_kds
geoms = [(0.05, 0.6), (0.001, 0.95), (0.05, 0.9877), (0.01, 0.99)]
for (Lam, a) in geoms:
    print(f'Lam={Lam} a={a}')
    for wf in (0.001, 0.01, 0.5, 0.99, 0.999):
        line = []
        for ell in (1.0, 10.0, 1e3, 1e5):
            r = check(Lam, a, wf, ell)
            if not r['signs']: line.append('SIGNS'); continue
            line.append(f"Phi={r['event_Phi_rel']:+.3f} C+={r['cosmo_C']:+.3f} L2'={r['cosmo_L2p']:+.3f}")
        print(f'   wfrac={wf}: ' + ' | '.join(line))

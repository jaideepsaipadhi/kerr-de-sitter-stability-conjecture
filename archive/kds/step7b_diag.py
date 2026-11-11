"""Diagnose the failures: which side fails, by how much, and where in the window."""
import numpy as np
from step7_verify import run_case
for (Lam, a) in [(0.05, 0.6), (0.02, 0.9), (0.1, 0.3), (0.001, 0.95)]:
    print(f'\n### Lam={Lam} a={a}')
    for m in (1, 2, 4):
        for lam in (2.0, 10.0, 50.0):
            row = []
            for frac in (0.05, 0.35, 0.65, 0.95):
                r = run_case(Lam, a, 1.0, m, lam, frac)
                if r['status'] == 'CLOSED':
                    row.append('  ok  '); continue
                tag = ''
                for side in ('cosmo', 'event'):
                    v = r.get(side, {})
                    if v.get('tau') is None or v['tau'] < 0: tag += side[0].upper() + 'lp '
                    elif 'maxC' in v:
                        if v['maxC'] > 1e-9: tag += side[0] + 'C '
                        if v['maxD'] > 1e-9: tag += side[0] + 'D '
                        if v['minf'] < -1e-9: tag += side[0] + 'f '
                if r.get('Lq_t', 1) <= 0: tag += 'Lq '
                row.append(tag.strip() or 'seam')
            print(f'  m={m} lam={lam:5.1f}: ' + ' | '.join(f'{x:>8}' for x in row))

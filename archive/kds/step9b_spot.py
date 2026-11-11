"""Spot checks of step 9 on geometries and modes it has not been tried on."""
from step9_freeq import run
cases = [(0.05, 0.9877, 1, 2.0, 0.05), (0.05, 0.9877, 2, 10.0, 0.5), (0.05, 0.9877, 4, 50.0, 0.95),
         (0.1, 0.3, 1, 2.0, 0.05), (0.1, 0.3, 4, 2.0, 0.95), (0.02, 0.9, 8, 120.0, 0.5),
         (0.08, 0.5, 2, 4.0, 0.5), (0.01, 0.9, 3, 4.0, 0.2)]
c = 0
for case in cases:
    st, info = run(*case); c += st == 'CLOSED'
    print(st, case, f"min-tau={info.get('best'):+.2e}", flush=True)
print(f'{c}/{len(cases)} closed')

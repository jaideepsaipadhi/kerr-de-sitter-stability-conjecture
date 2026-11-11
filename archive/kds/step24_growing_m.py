"""Step 24: does the free-qt method (step 9) keep closing as m grows, with lambda = m^2 ell?  Track the margin."""
import warnings; warnings.filterwarnings('ignore')
import sys, numpy as np
from step9_freeq import run
cases = [(0.05, 0.6, 0.5, 1.2), (0.001, 0.95, 0.2, 1.0), (0.01, 0.99, 0.5, 0.9), (0.05, 0.9877, 0.9, 0.8)]
ms = [int(x) for x in sys.argv[1:]] or [4, 16, 64]
for (Lam, a, wf, ell) in cases:
    line = []
    for m in ms:
        st, info = run(Lam, a, m, ell*m*m, wf)
        line.append(f'm={m}: {st} tau={info.get("best", float("nan")):+.1e}')
    print(f'Lam={Lam} a={a} wfrac={wf} ell={ell}: ' + ' | '.join(line), flush=True)

"""Step 26: sweep the explicit layer construction (step 25) over geometries, window positions, ell, and m (small m included)."""
import warnings; warnings.filterwarnings('ignore')
import io, contextlib, sys, numpy as np
with contextlib.redirect_stdout(io.StringIO()):
    from step17_ellmin import ell_min
from step25_layer import build, check
from step5_window import ctdc_window
from multiprocessing import Pool

def job(c):
    Lam, a, wf, ellmode, m = c
    lo, hi, _, _ = ctdc_window(Lam, a, 1.0)
    w = lo + wf*(hi - lo)
    ell = ell_min(Lam, a, w) if ellmode == 'min' else ellmode
    best = None
    for (delta, p) in ((0.5, 0.5), (1.0, 0.3), (0.25, 0.7), (2.0, 0.2), (0.5, 0.1)):
        try:
            r = check(build(Lam, a, m, ell, wf), delta, p)
        except Exception as ex:
            return c, False, str(ex)
        if r['ok']:
            return c, True, (delta, p)
        best = r
    return c, False, {k: round(v, 3) for k, v in best.items() if isinstance(v, float)}

if __name__ == '__main__':
    ms = [int(x) for x in sys.argv[1:]] or [1, 2, 4, 16, 64]
    geoms = [(0.05, 0.6), (0.02, 0.9), (0.1, 0.3), (0.001, 0.95), (0.05, 0.9877), (0.01, 0.99), (0.001, 0.999)]
    cases = [(L, a, wf, em, m) for (L, a) in geoms for wf in (0.02, 0.3, 0.7, 0.98) for em in ('min', 5.0) for m in ms]
    tally = {m: [0, 0] for m in ms}; fails = []
    with Pool() as pool:
        for c, ok, info in pool.imap_unordered(job, cases):
            tally[c[-1]][0] += ok; tally[c[-1]][1] += 1
            if not ok: fails.append((c, info))
    for m in ms: print(f'm={m}: {tally[m][0]}/{tally[m][1]} closed by the explicit layer construction')
    for f in sorted(fails, key=lambda x: x[0][-1])[:25]: print('  fail', f)

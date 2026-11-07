"""
Parallel certification driver for the event lemma on a compact core region.

Tiles: t in geometric bands, zf in bands, A in decades, Bf in bands.  Each tile goes to boxes2.certify_region
(adaptive splitting).  Results are appended to cert_log.txt as they arrive, so a run can be stopped and resumed:
tiles already marked OK are skipped.

Usage:
  python3 driver.py core      # t in [1.05, 300], zf in [0.02, 0.98], A in [1e-3, 1e6], Bf in [1, 100]
  python3 driver.py small     # quick sub-range, for testing
"""
import sys, os, time, itertools, json
from multiprocessing import Pool
from boxes2 import certify_region

LOG = 'cert_log.txt'


def tiles(spec):
    if spec == 'small':
        T = [1.2, 1.6, 2.5, 4.0]; Z = [0.05, 0.3, 0.7, 0.95]; A = [1e-2, 1, 1e2]; B = [1.0, 3.0, 10.0]
    else:
        T = [1.05, 1.2, 1.5, 2.0, 3.0, 5.0, 8.0, 13.0, 20.0, 35.0, 60.0, 100.0, 170.0, 300.0]
        Z = [0.02, 0.1, 0.3, 0.5, 0.7, 0.9, 0.98]
        A = [1e-3, 1e-2, 1e-1, 1, 10, 1e2, 1e3, 1e4, 1e5, 1e6]
        B = [1.0, 1.5, 3.0, 10.0, 100.0]
    for (t0, t1), (z0, z1), (a0, a1), (b0, b1) in itertools.product(zip(T, T[1:]), zip(Z, Z[1:]), zip(A, A[1:]), zip(B, B[1:])):
        yield ((t0, t1), (z0, z1), (a0, a1), (b0, b1))


def job(tile):
    t0 = time.time()
    done, fails = certify_region(tile, max_depth=16)
    return tile, done, fails, time.time() - t0


if __name__ == '__main__':
    spec = sys.argv[1] if len(sys.argv) > 1 else 'small'
    seen = set()
    if os.path.exists(LOG):
        for line in open(LOG):
            rec = json.loads(line)
            if rec['ok']: seen.add(json.dumps(rec['tile']))
    todo = [tl for tl in tiles(spec) if json.dumps(tl) not in seen]
    print(f'{len(todo)} tiles to do ({len(seen)} already certified)', flush=True)
    with Pool() as pool, open(LOG, 'a') as log:
        nok = nbad = 0
        for tile, done, fails, dt in pool.imap_unordered(job, todo):
            ok = len(fails) == 0
            nok += ok; nbad += not ok
            log.write(json.dumps(dict(tile=tile, ok=ok, boxes=done, fails=fails[:5], secs=round(dt, 1))) + '\n'); log.flush()
            print(f'{"OK " if ok else "BAD"} {tile}  boxes={done} fails={len(fails)} ({dt:.0f}s)   [{nok} ok, {nbad} bad]', flush=True)
    print(f'finished: {nok} tiles certified, {nbad} with uncertified boxes')

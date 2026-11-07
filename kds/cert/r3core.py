import itertools, sys, time
from multiprocessing import Pool
from r3cert import certify
T = [(1.01, 1.2), (1.2, 1.6), (1.6, 2.5), (2.5, 4), (4, 8), (8, 20), (20, 50), (50, 100)]
F = [(0.02, 0.2), (0.2, 0.5), (0.5, 0.8), (0.8, 0.98)]
A = [(1e-3, 1.0), (1.0, 1e2), (1e2, 1e4)]
Bf = [(1.0, 10.0), (10.0, 1e3)]
tiles = list(itertools.product(T, F, A, Bf))
def run(tile):
    d, f, n = certify(tile, tlimit=float(sys.argv[1]))
    return tile, d, len(f), n, f[:1]
if __name__ == '__main__':
    t0 = time.time(); tot = {'R1': 0, 'R2': 0, 'R3': 0}; bad = []
    with Pool(1) as P:
        for tile, d, nf, n, f in P.imap_unordered(run, tiles):
            for k in tot: tot[k] += d[k]
            if nf: bad.append((tile, nf, f))
            print(tile, d, 'uncert', nf, flush=True)
    print('TOTAL', len(tiles), 'tiles', tot, 'tiles with uncertified boxes:', len(bad), f'{time.time()-t0:.0f}s')
    for b in bad: print('  BAD', b)

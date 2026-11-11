"""R3 certification: interval-certify the event lemma with the ONE fixed multiplier (delta, beta, c) = (1, 2, 1/2),
skipping boxes covered analytically by R1/R2 (LEMMAS 47): a box is covered if, on all of it,
  2 zf^2 (t-1) (A/z0 + 3) < 1   (R2 condition; its points with p <= 1 are R1)   or   p = 2 zf (A + B - 3) <= 1  (R1)."""
import sys, time
import boxes2
from boxes2 import cert_box, Bfloor, _fwd, _inv
FIXED = (1.0, 2.0, 0.5)

def covered(box):
    (t0, t1), (f0, f1), (A0, A1), (b0, b1) = box
    z0lo = 1 + f0*(t0 - 1)
    if 2*f1**2*(t1 - 1)*(A1/z0lo + 3) < 1: return 'R2'
    if 2*f1*(A1 + b1*Bfloor(t1) - 3) <= 1: return 'R1'
    return None

def certify(region, max_depth=22, max_boxes=20000, tlimit=None):
    base = [max(_fwd(i, hi) - _fwd(i, lo), 1e-12) for i, (lo, hi) in enumerate(region)]
    todo = [(region, 0)]; done = {'R1': 0, 'R2': 0, 'R3': 0}; fails = []; tried = 0; T0 = time.time()
    while todo:
        box, d = todo.pop(); tried += 1
        if tried > max_boxes or (tlimit and time.time() - T0 > tlimit):
            fails.append(('CAP', box)); fails += [bx for bx, _ in todo]; break
        cv = covered(box)
        if cv: done[cv] += 1; continue
        if cert_box(box, FIXED): done['R3'] += 1; continue
        if d >= max_depth: fails.append(box); continue
        w = [(_fwd(i, hi) - _fwd(i, lo))/base[i] for i, (lo, hi) in enumerate(box)]
        i = max(range(4), key=lambda j: w[j]); lo, hi = box[i]
        mid = _inv(i, (_fwd(i, lo) + _fwd(i, hi))/2)
        b1 = list(box); b2 = list(box); b1[i] = (lo, mid); b2[i] = (mid, hi)
        todo += [(tuple(b1), d + 1), (tuple(b2), d + 1)]
    return done, fails, tried

if __name__ == '__main__':
    region = eval(sys.argv[1]); tl = float(sys.argv[2]) if len(sys.argv) > 2 else None
    t0 = time.time(); done, fails, tried = certify(region, tlimit=tl)
    print(region, 'certified', done, 'uncertified', len(fails), 'boxes tried', tried, f'{time.time()-t0:.0f}s', fails[:2], flush=True)

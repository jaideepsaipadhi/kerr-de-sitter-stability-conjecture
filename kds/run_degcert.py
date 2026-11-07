"""Certify EXACT deg_sigma E_i = 29 at all three test points (closes the last gap of the (sigma, m) argument, LEMMAS 69).
For each point and minor: exact K=9 reversed-Toeplitz degree bound at mu0 = 1..(D_mu + 1), D_mu = proved mu-degree bound
(C 61, G1 7, G2 10).  If max over mu0 equals the per-minor target (C 11, G1 6, G2 12; sum 29 = mod-p degree), every
sigma-coefficient above the target vanishes at D_mu + 1 points, hence identically.  2 workers; ~50 min on 2 cores."""
import subprocess, re, sys
from concurrent.futures import ThreadPoolExecutor
D = {'C': 62, 'G1': 8, 'G2': 11}; T = {'C': 11, 'G1': 6, 'G2': 12}
jobs = [(ip, nm) for nm in ('C', 'G1', 'G2') for ip in (0, 1, 2)]
def run(job):
    ip, nm = job
    mus = ','.join(str(m) for m in range(1, D[nm] + 1))
    out = subprocess.run([sys.executable, 'step58i_fixed_mu_deg.py', str(ip), '9', mus, nm], capture_output=True, text=True).stdout
    degs = [int(x) for x in re.findall(r'deg_sigma <= (\d+)', out)]
    open(f'degcert_pt{ip}_{nm}.log', 'w').write(out)
    return ip, nm, len(degs), max(degs) if degs else None
with ThreadPoolExecutor(2) as ex:
    res = list(ex.map(run, jobs))
ok = True
for ip, nm, n, mx in sorted(res):
    good = (n == D[nm]) and (mx is not None) and (mx <= T[nm]); ok &= good
    print(f'pt{ip} {nm}: {n}/{D[nm]} mu0 values, max exact deg bound {mx} (target {T[nm]})  {"OK" if good else "FAIL"}')
print('EXACT deg_sigma E_i <= 29 at all three points:', ok)

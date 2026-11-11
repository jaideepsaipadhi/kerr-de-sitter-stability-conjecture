import numpy as np, time
exec(open('phys149.py').read().split("rng=np.random.default_rng(1)")[0])
rng=np.random.default_rng(8); t0=time.time(); up=[]; lo=[]
while time.time()-t0<45:
    r0=1-10**rng.uniform(-3,0)*0.999; r2=10**rng.uniform(0.05,3); frac=rng.uniform(0.001,0.99999); m=int(rng.integers(1,6))
    try: tt,zz,AA,BB=pdata(r0,1.0,r2,m,frac)
    except Exception: continue
    x0=zz-1; zs=1+0.62*x0*x0/zz
    q=1+(zs-1)*np.linspace(1e-4,1,200); s=zz-q; L=AA/q+BB-3*q; Lp=AA/q**2+3; mu=q*(q-1)*(tt-q)
    g=s*s*L*L+mu*(L-s*Lp); N=Nf(tt,zz,AA,BB,s)
    # normalise: kappa in units of g-scale/N-scale: use scale ratio  (s^2 L^2 + mu L)/|N| ... record r=g/N
    sc=(s*s*L*L+mu*L)/np.abs(N)
    r=g/N/sc
    if np.any(N>0): up.append(r[N>0].min())
    m_=(N<0)&(g<=0)
    if np.any(m_): lo.append(r[m_].max())
print('need kappa_n <',min(up) if up else None,' and >',max(lo) if lo else None, len(up),len(lo))

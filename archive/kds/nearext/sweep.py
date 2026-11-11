import numpy as np, time
exec(open('/tmp/ev2.py').read().split("if __name__")[0])
def cosmo(d,kap,n=3000):
    t,A,B,z0,q1=d['t'],d['A'],d['B'],d['z0'],abs(d['q1'])
    mu=lambda z: z*(z-1)*(t-z); mup=lambda z: (z-1)*(t-z)+z*(t-z)-z*(z-1)
    q=lambda z: q1*(z-z0)+kap*(z-1)*(t-z); qp=lambda z: q1+kap*(t+1-2*z)
    L=lambda z: A/z+B-3*z+(kap/z)*(2*q1*(z-z0)+kap*(z-1)*(t-z))
    zt=brentq(q,1+1e-12,t-1e-9); zs=zt+(t-zt)*np.linspace(1e-4,1-1e-6,n); e=1e-7
    Lp=(L(zs*(1+e))-L(zs*(1-e)))/(2*e*zs)
    return bool((mup(zs)*q(zs)-2*mu(zs)*qp(zs)).max()<=0 and (q(zs)*Lp-kap*mup(zs)).max()<=0 and L(t)>0 and qp(zs).min()>0)
rng=np.random.default_rng(1); t0=time.time(); bad=[]; N=0
while time.time()-t0<50:
    r0=1-10**rng.uniform(-2.5,0); r2=10**rng.uniform(0.1,2.5); frac=rng.uniform(0.02,0.999); m=int(rng.integers(1,4))
    try: d=data(r0,1.0,r2,m,frac)
    except Exception: continue
    k0=-d['Lt']*d['t']/(2*d['q1']*(d['t']-d['z0'])); kap=max(0.0,1.5*k0)
    try:
        ev=event(d,kap,n=2000); co=cosmo(d,kap)
    except Exception as ex: bad.append(('err',r0,r2,frac,m,str(ex)[:40])); continue
    N+=1
    ok= co and isinstance(ev,dict) and ev['changes']==1 and ev['neg_start'] and ev['fpos'] and ev['seam']
    if not ok: bad.append((round(r0,4),round(r2,2),round(frac,3),m,round(kap,3),co,ev))
print('cases',N,'bad',len(bad)); [print(b) for b in bad[:12]]

import numpy as np, time
exec(open('/tmp/chk.py').read().split("for r0 in (0.9")[0])
exec(open('/tmp/evab.py').read().split("rng=")[0])
rng=np.random.default_rng(11); t0=time.time(); rows=[]
while time.time()-t0<40:
    r0=1-10**rng.uniform(-2.5,0); r2=10**rng.uniform(0.1,2.5); frac=rng.uniform(0.02,0.999); m=int(rng.integers(1,6))
    try: d=data(r0,1.0,r2,m,frac)
    except Exception: continue
    t,z0,A,B=d['t'],d['z0'],d['A'],d['B']; x0=z0-1
    Lz0=A/z0+B-3*z0; p=2*x0*(A+B-3)/(t-1)
    zz=1+x0*np.linspace(1e-6,1-1e-6,3000); Nv=Nf(zz,t,z0,A,B); ch=np.sum(np.diff(np.sign(Nv))!=0)
    K=None
    if ch==1:
        i=np.where(np.diff(np.sign(Nv))!=0)[0][0]; zh=zz[i]; sh=z0-zh
        H=2*zh*(zh-1)*(t-zh)+sh*((zh-1)*(t-zh)+zh*(t-zh)-zh*(zh-1)); K=H-sh**3*(A/zh**2+3)
    rows.append((Lz0/(abs(A)+abs(B)+1),p,ch,K, B-(5.5*t-1.25), d['lam']))
r=np.array([(a,b,c,(k if k is not None else np.nan),e,f) for a,b,c,k,e,f in rows])
print('N',len(r),'min Lz0 rel',r[:,0].min(),'#ch!=1 (p>1)',np.sum((r[:,2]!=1)&(r[:,1]>1)),'#K<=0 (p>1)',np.sum((r[:,3]<=0)&(r[:,1]>1)),'#below floor',np.sum(r[:,4]<0),'min lam',r[:,5].min())

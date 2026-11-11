exec(open('step138.py').read().split("rng=np.random")[0])
from scipy.integrate import cumulative_trapezoid as ct
rng=np.random.default_rng(4); tot=fail=0; mins=[]
for i in range(3000):
    tt=1+10**rng.uniform(-2,3); zz=1+(tt-1)*rng.uniform(1e-2,1-1e-2); x0=zz-1
    a0=max(0.0,zz*((tt-1)/(2*x0**2)-3))
    AA=a0+10**rng.uniform(-3,4)*rng.integers(0,2)*(a0+1)
    bb=CC*AA/zz+DD*zz+10**rng.uniform(-4,1); BB=3*zz-AA/zz+bb
    if 2*x0*(AA+BB-3)<=tt-1: continue
    c=[mp.mpf(g(tt,zz,AA,BB)) for g in cf]
    r=[mp.re(x) for x in mp.polyroots(c,maxsteps=300,extraprec=300) if abs(mp.im(x))<1e-15 and 0<mp.re(x)<x0]
    Sh=float(r[0]); zh=zz-Sh; mu_=zh*(zh-1)*(tt-zh); mup=(zh-1)*(tt-zh)+zh*(tt-zh)-zh*(zh-1)
    H=2*mu_+Sh*mup; w0=H/Sh
    q=zh+Sh*(1-np.concatenate([np.geomspace(1,1e-9,40000),[0]]))[::-1][::-1]; q=np.sort(zh+Sh*np.concatenate([[0],np.geomspace(1e-12,1,60000)])); sq=zz-q; Lp=AA/q**2+3; muq=q*(q-1)*(tt-q); Lq=AA/q+BB-3*q
    ok=False
    for fac in (1,1.1,1.25,1.5,2,3):
        M=fac*w0; I=ct((2*sq*Lq-M)/muq,q,initial=0); v=np.exp(I)*(w0-ct(2*sq*Lp*np.exp(-I),q,initial=0))
        if v.max()<=M*(1+1e-9) and v[-1]>0: ok=True; best=v[-1]/w0; break
    tot+=1; fail+= not ok
    if not ok: print(round(tt,3),round(x0/(tt-1),3),AA,bb,round(Sh/x0,3))
print('tot',tot,'fail',fail)

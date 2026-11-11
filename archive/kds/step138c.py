exec(open('step138.py').read().split("rng=np.random")[0])
rng=np.random.default_rng(4); res={a:[0,0] for a in (0.02,0.05,0.1,0.2,0.3)}
for i in range(3000):
    tt=1+10**rng.uniform(-2,3); zz=1+(tt-1)*rng.uniform(1e-2,1-1e-2); x0=zz-1
    a0=max(0.0,zz*((tt-1)/(2*x0**2)-3))
    AA=a0+10**rng.uniform(-3,4)*rng.integers(0,2)*(a0+1)
    bb=CC*AA/zz+10**rng.uniform(-4,1); BB=3*zz-AA/zz+bb
    if 2*x0*(AA+BB-3)<=tt-1: continue
    c=[mp.mpf(g(tt,zz,AA,BB)) for g in cf]
    r=[mp.re(x) for x in mp.polyroots(c,maxsteps=300,extraprec=300) if abs(mp.im(x))<1e-15 and 0<mp.re(x)<x0]
    Sh=float(r[0]); zh=zz-Sh; mu_=zh*(zh-1)*(tt-zh); mup=(zh-1)*(tt-zh)+zh*(tt-zh)-zh*(zh-1)
    H=2*mu_+Sh*mup; w0=H/Sh
    q=np.linspace(zh,zz,3001); sq=zz-q; Lq=AA/q+BB-3*q; Lp=AA/q**2+3; muq=q*(q-1)*(tt-q)
    for a in res:
        k=w0*(a+(1-a)*sq/Sh)
        cond=k**2+2*muq*sq*Lp-2*sq*Lq*k-muq*(1-a)*w0/Sh
        res[a][0]+=1; res[a][1]+= np.max(cond)>0
print(res)

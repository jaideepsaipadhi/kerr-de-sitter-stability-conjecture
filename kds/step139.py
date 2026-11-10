exec(open('step138.py').read().split("rng=np.random")[0])
rng=np.random.default_rng(4); fams=[(0.0,2.0)]
anyok=0; tot=0; cnt={f:0 for f in fams}
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
    sg=np.linspace(0,1,2001); q=zz-Sh*sg; sq=Sh*sg; Lq=AA/q+BB-3*q; Lp=AA/q**2+3; muq=q*(q-1)*(tt-q)
    tot+=1; ok=False
    for (cc,p) in fams:
        k=w0*(cc+(1-cc)*sg**p); kz=w0*(1-cc)*p*sg**(p-1)*(-1/Sh) if p>=1 else w0*(1-cc)*p*np.where(sg>0,sg,1e-12)**(p-1)*(-1/Sh)
        cond=muq*kz+k**2-2*sq*Lq*k+2*muq*sq*Lp
        if np.max(cond)<=0: cnt[(cc,p)]+=1; ok=True
    anyok+=ok
print(tot,anyok,sorted(cnt.items(),key=lambda x:-x[1])[:5])

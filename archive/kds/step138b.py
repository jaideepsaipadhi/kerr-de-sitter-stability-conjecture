exec(open('step138.py').read().split("rng=np.random")[0])
from scipy.integrate import solve_ivp
rng=np.random.default_rng(4); k=0
for i in range(3000):
    tt=1+10**rng.uniform(-2,3); zz=1+(tt-1)*rng.uniform(1e-2,1-1e-2); x0=zz-1
    a0=max(0.0,zz*((tt-1)/(2*x0**2)-3))
    AA=a0+10**rng.uniform(-3,4)*rng.integers(0,2)*(a0+1)
    bb=CC*AA/zz+10**rng.uniform(-4,1); BB=3*zz-AA/zz+bb
    if 2*x0*(AA+BB-3)<=tt-1: continue
    c=[mp.mpf(g(tt,zz,AA,BB)) for g in cf]
    r=[mp.re(x) for x in mp.polyroots(c,maxsteps=300,extraprec=300) if abs(mp.im(x))<1e-15 and 0<mp.re(x)<x0]
    Sh=float(r[0]); zh=zz-Sh; mu_=zh*(zh-1)*(tt-zh); mup=(zh-1)*(tt-zh)+zh*(tt-zh)-zh*(zh-1)
    H=2*mu_+Sh*mup; Lf=lambda q: AA/q+BB-3*q
    def rhs(q,w):
        sq=zz-q; muq=q*(q-1)*(tt-q); return [(-w[0]**2+2*sq*Lf(q)*w[0]-2*muq*sq*(AA/q**2+3))/muq]
    sol=solve_ivp(rhs,(zh,zz),[H/Sh],rtol=1e-10,atol=1e-12,dense_output=True)
    pts=[0,0.25,0.5,0.75,0.9,1.0]
    print(round(tt,2),round(x0/(tt-1),2),'Sh/x0',round(Sh/x0,3),'w/w0 along:',[round(float(sol.sol(zh+p*Sh)[0])/(H/Sh),3) for p in pts], 'L0*mu0', round(bb*zz*(zz-1)*(tt-zz)/(H/Sh)**2,3))
    k+=1
    if k>10: break

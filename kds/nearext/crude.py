import numpy as np, time
exec(open('/tmp/chk.py').read().split("for r0 in (0.9")[0])
def crude(m,al,xi):
    l0=max(abs(m),2); return l0*(l0+1)-4-4*al-4*al*abs(m)*(1+al)-2*(1+al)**2*abs(xi)*(abs(m)+2)
rng=np.random.default_rng(13); t0=time.time(); worst=[]; n=0
while time.time()-t0<40:
    r0=1-10**rng.uniform(-3,0); r2=10**rng.uniform(0.05,3); frac=rng.uniform(0,1); m=int(rng.integers(1,12))
    try: d=data(r0,1.0,r2,m,frac)
    except Exception: continue
    # recompute B with crude lambda
    t,z0,A=d['t'],d['z0'],d['A']; C0=d['C0']; c1=d['B']-d['lam']*C0
    L_,a,al,w1=kds(r0,1.0,r2); r3=-(r0+1+r2)
    w2=None
    om=None
    # xi = a*sigma: recover from data? recompute window w
    Dp=lambda r: -(L_/3)*sum(np.prod([r-x for k,x in enumerate([r0,1,r2,r3]) if k!=j]) for j in range(4))
    a2=a*a; P=[(x*x+a2)/abs(Dp(x)) for x in [r0,1,r2,r3]]; Dabs=[abs(Dp(x)) for x in [r0,1,r2,r3]]
    w2=a*(1/Dabs[1]-1/Dabs[0])/(P[1]-P[0]); w=w2+frac*(w1-w2); xi=a*m*w
    Bc=crude(m,al,xi)*C0+c1; Lz0=A/z0+Bc-3*z0
    n+=1; worst.append((Lz0/(A/z0+abs(Bc)+3*z0),round(r0,4),round(r2,2),round(frac,3),m))
worst.sort(); print(n); [print(w) for w in worst[:6]]

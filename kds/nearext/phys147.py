import numpy as np, time, sympy as sp
exec(open('/tmp/chk.py').read().split("for r0 in (0.9")[0])
z,t,z0,A,B,S=sp.symbols('z t z0 A B S')
s_=z0-z; mus=z*(z-1)*(t-z); Ls=A/z+B-3*z; F=mus/s_**2
nu,dn=sp.fraction(sp.factor(sp.together(sp.diff(mus*sp.diff(F,z),z)-2*s_*sp.diff(Ls*F,z))))
Nf=sp.lambdify((t,z0,A,B,S),sp.expand(nu.subs(z,z0-S)))
rng=np.random.default_rng(1); t0=time.time(); stats={'R3':0,'ok':0,'tot':0}; fails=[]
ths=np.linspace(0.3,1.95,34)
while time.time()-t0<50:
    r0=1-10**rng.uniform(-3,0)*0.999; r2=10**rng.uniform(0.05,3); frac=rng.uniform(0.001,0.999); m=int(rng.integers(1,30))
    try: d=data(r0,1.0,r2,m,frac)
    except Exception: continue
    tt,zz,AA,BB=d['t'],d['z0'],d['A'],d['B']; x0=zz-1; stats['tot']+=1
    if not (2*x0*(AA+BB-3)>tt-1 and 2*x0**2*(AA/zz+3)>=tt-1): continue
    stats['R3']+=1
    u=np.linspace(1e-6,1-1e-6,800); q=1+x0*u; s=zz-q; L=AA/q+BB-3*q; Lp=AA/q**2+3; mu=q*(q-1)*(tt-q)
    mup=(q-1)*(tt-q)+q*(tt-q)-q*(q-1); H=2*mu+s*mup; N=Nf(tt,zz,AA,BB,s)
    good=False
    for th in ths:
        g=th*(2-th)*s*s*L*L+th*mu*L-(2-th)*mu*s*Lp
        if not np.any((N>=0)&((g<=0)|(H-th*s*s*L<0))): good=True; break
    stats['ok']+=good
    if not good: fails.append((r0,r2,frac,m))
print(stats, fails[:5])

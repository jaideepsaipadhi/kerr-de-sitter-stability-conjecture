import numpy as np, time, sympy as sp
exec(open('/tmp/win.py').read().split("worst=[]")[0])
z,t,z0,A,B,S=sp.symbols('z t z0 A B S')
s_=z0-z; mus=z*(z-1)*(t-z); Ls=A/z+B-3*z; F=mus/s_**2
nu,dn=sp.fraction(sp.factor(sp.together(sp.diff(mus*sp.diff(F,z),z)-2*s_*sp.diff(Ls*F,z))))
Nf=sp.lambdify((t,z0,A,B,S),sp.expand(nu.subs(z,z0-S)))
def pdata(r0,r1,r2,m,frac):
    r3=-(r0+r1+r2); rs=[r0,r1,r2,r3]
    L,a,al,w1=kds(r0,r1,r2); Xi=1+al; a2=a*a
    Dp=lambda r: -(L/3)*sum(np.prod([r-x for k,x in enumerate(rs) if k!=j]) for j in range(4))
    P=[(x*x+a2)/abs(Dp(x)) for x in rs]; Dabs=[abs(Dp(x)) for x in rs]
    w2=a*(1/Dabs[1]-1/Dabs[0])/(P[1]-P[0]); w=w2+frac*(w1-w2); om=m*w
    kap=[1j*Xi*(om*(x*x+a2)-a*m)/Dp(x) for x in rs]
    t=(r1-r3)*(r2-r0)/((r1-r0)*(r2-r3)); X=[abs(kap[j]+kap[3]) for j in (0,1,2)]
    A=t*(X[0]**2-(X[1]+X[2])**2); C0=3/(L*(r1-r0)*(r2-r3)); c1=1+t+3*(1-al)*C0
    xi=a*om; l0=max(m,2)
    lamc=l0*(l0+1)-4-4*al-4*al*m*(1+al)-2*(1+al)**2*abs(xi)*(m+2)
    e1=2*(kap[1]+kap[3]); et=2*(kap[2]+kap[3]); q1v=(t-1)*e1.imag/2; qtv=-t*(t-1)*et.imag/2
    z0=1+(t-1)*q1v/(q1v-qtv)
    return t,z0,A,lamc*C0+c1
rng=np.random.default_rng(1); t0=time.time(); st={'tot':0,'R3':0,'Lz0neg':0,'a_ok':0,'ab_ok':0}; per={round(th,2):0 for th in np.linspace(0.3,1.95,34)}
ths=np.linspace(0.3,1.95,34)
while time.time()-t0<40:
    r0=1-10**rng.uniform(-3,0)*0.999; r2=10**rng.uniform(0.05,3); frac=rng.uniform(0.001,0.999); m=int(rng.integers(1,60))
    try: tt,zz,AA,BB=pdata(r0,1.0,r2,m,frac)
    except Exception: continue
    x0=zz-1; st['tot']+=1
    if AA/zz+BB-3*zz<=0: st['Lz0neg']+=1; continue
    if not (2*x0*(AA+BB-3)>tt-1 and 2*x0**2*(AA/zz+3)>=tt-1): continue
    st['R3']+=1
    u=np.linspace(1e-6,1-1e-6,800); q=1+x0*u; s=zz-q; L=AA/q+BB-3*q; Lp=AA/q**2+3; mu=q*(q-1)*(tt-q)
    mup=(q-1)*(tt-q)+q*(tt-q)-q*(q-1); H=2*mu+s*mup; N=Nf(tt,zz,AA,BB,s)
    ga=gb=False
    for th in ths:
        g=th*(2-th)*s*s*L*L+th*mu*L-(2-th)*mu*s*Lp
        a_=not np.any((N>=0)&(g<=0)); b_=not np.any((N>=0)&(H-th*s*s*L<0))
        ga|=a_
        if a_ and b_: gb=True; per[round(th,2)]+=1
    st['a_ok']+=ga; st['ab_ok']+=gb
print(st); print({k:v for k,v in per.items() if v})

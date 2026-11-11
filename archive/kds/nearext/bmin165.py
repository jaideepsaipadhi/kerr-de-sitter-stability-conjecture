import numpy as np, time, sympy as sp, warnings; warnings.filterwarnings('ignore')
exec(open('phys149.py').read().split("rng=np.random.default_rng(1)")[0])
z_,t_,z0_,A_,B_=sp.symbols('z t z0 A B')
s=z0_-z_; mu=z_*(z_-1)*(t_-z_); L=A_/z_+B_-3*z_; F=mu/s**2
nu,dn=sp.fraction(sp.factor(sp.together(sp.diff(mu*sp.diff(F,z_),z_)-2*s*sp.diff(L*F,z_))))
N=sp.expand(nu); H=2*mu+s*sp.diff(mu,z_); T=sp.expand(H-s**2*L)
NB=N.coeff(B_); Ts=sp.together(T+s**2*N/NB); n,d=sp.fraction(sp.factor(Ts)); n=sp.expand(n)
AT=sp.lambdify((t_,z0_,z_),-n.subs(A_,0)/n.coeff(A_))
Nf=sp.lambdify((t_,z0_,A_,B_,z_),N)
def pdata2(r0,r1,r2,m,frac):
    r3=-(r0+r1+r2); rs=[r0,r1,r2,r3]
    L,a,al,w1=kds(r0,r1,r2); Xi=1+al; a2=a*a
    Dp=lambda r: -(L/3)*sum(np.prod([r-x for k,x in enumerate(rs) if k!=j]) for j in range(4))
    P=[(x*x+a2)/abs(Dp(x)) for x in rs]; Dabs=[abs(Dp(x)) for x in rs]
    w2=a*(1/Dabs[1]-1/Dabs[0])/(P[1]-P[0]); w=w2+frac*(w1-w2); om=m*w
    kap=[1j*Xi*(om*(x*x+a2)-a*m)/Dp(x) for x in rs]
    t=(r1-r3)*(r2-r0)/((r1-r0)*(r2-r3)); X=[abs(kap[j]+kap[3]) for j in (0,1,2)]
    A=t*(X[0]**2-(X[1]+X[2])**2); C0=3/(L*(r1-r0)*(r2-r3)); c1=1+t+3*(1-al)*C0
    b2=C0*(1-2*Xi*Xi*a*w); b1=C0*(1-4*al*Xi-4*Xi*Xi*a*w); b0=C0*(-4-4*al)+c1
    Blow=m*m*(b2+min(b1,0))+b0+(max(b1,0)*m if False else 0)
    if m==1: Blow=b2+b1+b0+4*C0
    e1=2*(kap[1]+kap[3]); et=2*(kap[2]+kap[3]); q1v=(t-1)*e1.imag/2; qtv=-t*(t-1)*et.imag/2
    z0=1+(t-1)*q1v/(q1v-qtv)
    return t,z0,A,Blow,b1
NEGB1=[0]
rng=np.random.default_rng(3); t0=time.time(); R=[]
while time.time()-t0<40:
    r0=1-10**rng.uniform(-3.5,0)*0.999; r2=10**rng.uniform(0.02,3.5); frac=rng.uniform(0.0001,0.99999); m=int(rng.integers(1,40))
    try: tt,zz,AA,BB,b1=pdata2(r0,1.0,r2,m,frac); NEGB1[0]+=b1<0
    except Exception: continue
    x0=zz-1
    if not (2*x0*(AA+BB-3)>tt-1 and 2*x0**2*(AA/zz+3)>=tt-1): continue
    q=1+x0*np.linspace(1e-5,1-1e-3,800); at=AT(tt,zz,q)
    i=np.where(at>=AA)[0]
    if len(i)==0: R.append((np.inf,x0,tt)); continue
    zT=q[i.min()]; Nv=Nf(tt,zz,AA,BB,zT)
    R.append((Nv,x0,tt))
R=np.array(R); print('b1<0 count',NEGB1[0]); print(len(R),'N(zT;Bc)>0 (violations):',(R[:,0]>0).sum(), ' zT exists frac',np.isfinite(R[:,0]).mean())

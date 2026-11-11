import sympy as sp, numpy as np
z,t,z0,A,B,g=sp.symbols('z t z0 A B gamma')
s=z0-z; mu=z*(z-1)*(t-z); L=A/z+B-3*z; F=mu/s**2
nu,dn=sp.fraction(sp.factor(sp.together(sp.diff(mu*sp.diff(F,z),z)-2*s*sp.diff(L*F,z))))
N=sp.expand(nu); H=2*mu+s*sp.diff(mu,z); T=sp.expand(H-s**2*L)
Bsub=A/(2*z0)+2*t+g-A/z0+3*z0
Ng=sp.expand(N.subs(B,Bsub)); Tg=sp.expand(T.subs(B,Bsub))
f=sp.lambdify((t,z0,z),[Tg.coeff(g),Ng.coeff(g),Tg.coeff(A),Ng.coeff(A),Tg.subs({A:0,g:0}),Ng.subs({A:0,g:0})])
rng=np.random.default_rng(1); n=500000
x0=10**rng.uniform(-4,4,n); zz=1+x0; tt=zz+x0*10**rng.uniform(-5,6,n); q=1+x0*rng.uniform(1e-6,1-1e-6,n)
Tg_,Ng_,TA,NA,T0,N0=f(tt,zz,q)
kmin=Tg_/Ng_  # since cgamma = Tg - k Ng >=0 with Ng<0 -> k >= Tg/Ng
# cA = TA - k NA >=0 ; c0 = T0 - k N0 >=0
lo=np.maximum(kmin,0.0); hi=np.full(n,np.inf)
a0=np.maximum(0,zz*((tt-1)/(2*x0**2)-3))
# p-line: pA*A + pg*g >= rhs
pA=2*x0*(x0/zz+1/(2*zz)); pg=2*x0; rhs=(tt-1)-2*x0*(2*tt+3*x0)
verts=[]
# vertex1: A=a0, g=max(0,(rhs-pA*a0)/pg)
g1=np.maximum(0,(rhs-pA*a0)/pg); verts.append((a0,g1))
# vertex2: g=0, A=max(a0,rhs/pA)
A2=np.maximum(a0,rhs/pA); verts.append((A2,np.zeros(n)))
for Tc,Nc in [(TA,NA)]+[(T0+TA*va+Tg_*vg,N0+NA*va+Ng_*vg) for va,vg in verts]:
    # Tc - k Nc >= 0
    m=Nc>0; hi[m]=np.minimum(hi[m],Tc[m]/Nc[m])
    m=Nc<0; lo[m]=np.maximum(lo[m],Tc[m]/Nc[m])
    m=Nc==0; hi[m&(Tc<0)]=-np.inf
print('Ng<0 always',(Ng_<0).all(),' infeasible',(lo>hi).sum(),'of',n)
b=lo>hi
if b.any(): print('x0',np.percentile(x0[b],[0,50,100]),'t/z0',np.percentile((tt/zz)[b],[0,50,100]),'u',np.percentile(((q-1)/x0)[b],[0,50,100]))

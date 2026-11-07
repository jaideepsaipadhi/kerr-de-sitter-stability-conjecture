import sympy as sp, numpy as np, sys
from scipy.linalg import eig
exec(open('qnm_scalar2.py').read().split('import pickle')[0])
def solve(sv,av,Qv,mv,lamv,N):
    c=build(sv,av,Qv,mv,lamv); f=[sp.lambdify((u,w),sp.expand(ci*(1-u)**2),'numpy') for ci in c]
    k=np.arange(N); x=np.cos(np.pi*(2*k+1)/(2*N)); uu=(x+1)/2
    # Chebyshev basis T_n(2u-1) and derivatives at nodes
    T=np.polynomial.chebyshev
    V=np.array([T.chebval(x,np.eye(N)[n]) for n in range(N)]).T
    V1=np.array([T.chebval(x,T.chebder(np.eye(N)[n]))*2 for n in range(N)]).T
    V2=np.array([T.chebval(x,T.chebder(np.eye(N)[n],2))*4 for n in range(N)]).T
    def Mw(wv): return sum(np.diag(np.broadcast_to(f[j](uu,wv),uu.shape).astype(complex))@(V,V1,V2)[j] for j in range(3))
    M0,M1,Mm,M2=Mw(0),Mw(1),Mw(-1),Mw(2)
    A0=M0;A2=(M1+Mm-2*M0)/2;A1=(M1-Mm)/2
    print(" deg chk",np.abs(M2-(A0+2*A1+4*A2)).max()/np.abs(A0).max())
    n=N;Z=np.zeros((n,n));I=np.eye(n)
    ev=eig(np.block([[Z,I],[-A0,-A1]]),np.block([[I,Z],[Z,A2]]),right=False)
    ev=ev[np.isfinite(ev)&(abs(ev)<3)&(ev.real>0.05)]
    return sorted(ev,key=lambda z:-z.imag)
for N in (30,45,60):
    print(N,[f"{z:.6f}" for z in solve(-2,sp.Rational(1,10**6),0,2,4,N)[:3]])

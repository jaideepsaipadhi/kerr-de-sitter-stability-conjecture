# Spectral QNM solver (radial Teukolsky, spin s, on KN background with separation constant lam) -- validation step.
# R = P h, P = (r-r+)^(-s-i sig) (r-r-)^(-1-s+i sig+2iMw) e^{iwr}; h regular on u=1-r+/r in [0,1]; Chebyshev collocation;
# operator is polynomial (deg<=2) in w -> companion linearisation.
import numpy as np
from scipy.linalg import eig
def cheb(N):
    x=np.cos(np.pi*np.arange(N+1)/N); c=np.hstack([2,np.ones(N-1),2])*(-1)**np.arange(N+1)
    X=np.tile(x,(N+1,1)).T; dX=X-X.T; D=np.outer(c,1/c)/(dX+np.eye(N+1)); D-=np.diag(D.sum(1)); return D,x
def opmat(w,s,a,Q,m,lam,N,M=1.0):
    rp=M+np.sqrt(M*M-a*a-Q*Q); rm=M-np.sqrt(M*M-a*a-Q*Q)
    D,x=cheb(N); u=(x+1)/2; Du=2*D; Du2=Du@Du
    r=rp/(1-u+1e-300); r[-1]=1e12  # u=1 -> infinity handled by mapping; use large r (regularity)
    r=np.where(u>1-1e-14,1e8,r)
    OH=a/(rp*rp+a*a); sig=(w-m*OH)*(rp*rp+a*a)/(rp-rm)
    A=-s-1j*sig; B=-1-s+1j*sig+2j*M*w
    Dl=(r-rp)*(r-rm); K=(r*r+a*a)*w-a*m
    p1=A/(r-rp)+B/(r-rm)+1j*w; dp1=-A/(r-rp)**2-B/(r-rm)**2
    du=rp/r**2; d2u=-2*rp/r**3   # du/dr, d2u/dr2
    # h' = du h_u ; h'' = du^2 h_uu + d2u h_u
    c2=Dl**2; c1=(s+1)*Dl*(2*r-2*M); c0=K*K-2j*s*(r-M)*K+Dl*(4j*s*w*r-lam)
    # divide whole eq by r^2-scale to tame infinity: use weight 1/r^4
    wgt=1/(r**4)
    T2=(c2*du**2)*wgt; T1=(c2*(d2u+2*p1*du)+c1*du)*wgt; T0=(c2*(dp1+p1*p1)+c1*p1+c0)*wgt
    return np.diag(T2)@Du2+np.diag(T1)@Du+np.diag(T0)
def qnms(s,a,Q,m,lam,N=60):
    W=[0,1,-1,2]; Ms=[opmat(w,s,a,Q,m,lam,N) for w in W]
    A0=Ms[0]; A2=(Ms[1]+Ms[2]-2*Ms[0])/2; A1=(Ms[1]-Ms[2])/2
    print("deg check (should be ~0):",np.abs(Ms[3]-(A0+2*A1+4*A2)).max()/np.abs(A0).max())
    n=A0.shape[0]; Z=np.zeros((n,n)); I=np.eye(n)
    ev=eig(np.block([[Z,I],[-A0,-A1]]),np.block([[I,Z],[Z,A2]]),right=False)
    ev=ev[np.isfinite(ev)&(np.abs(ev)<3)&(ev.imag<0)&(ev.real>0.05)]
    return sorted(ev,key=lambda z:-z.imag)
for N in (40,60):
    e=qnms(-2,1e-6,0.0,2,4.0,N); print(N,[f"{z:.6f}" for z in e[:3]])

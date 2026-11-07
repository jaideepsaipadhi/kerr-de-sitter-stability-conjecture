# Geometry-only probe on extremal Kerr-Newman (M=1, a^2+Q^2=1): spin-weight -2 Teukolsky-type test field,
# Einstein-Maxwell COUPLING DROPPED. Same near/far matching as balance_kerr.py.
import numpy as np, sys
from scipy.integrate import solve_ivp
from sympy.physics.wigner import clebsch_gordan as CG
s=int(sys.argv[2]) if len(sys.argv)>2 else -2; m=int(sys.argv[1]) if len(sys.argv)>1 else 2
def specA(c):
    ls=list(range(max(abs(s),abs(m)),m+40)); n=len(ls); C=np.zeros((n,n))
    for i,l in enumerate(ls):
        for j,lp in enumerate(ls):
            if abs(l-lp)<=1: C[i,j]=float(np.sqrt((2*lp+1)/(2*l+1))*CG(lp,1,l,m,0,m)*CG(lp,1,l,-s,0,-s))
    H=np.diag([l*(l+1)-s*(s+1) for l in ls])-c*c*(C@C)+2*c*s*C
    return np.sort(np.linalg.eigvals(H).real)[0]
def integ(f,y0,x0,x1):
    def F(x,y):
        z=y[0]+1j*y[1]; w=y[2]+1j*y[3]; dz,dw=f(x,z,w); return [dz.real,dz.imag,dw.real,dw.imag]
    r=solve_ivp(F,(x0,x1),[y0[0].real,y0[0].imag,y0[1].real,y0[1].imag],rtol=1e-11,atol=1e-14,method='DOP853')
    return r.y[0,-1]+1j*r.y[1,-1], r.y[2,-1]+1j*r.y[3,-1]
for a in [float(v) for v in sys.argv[2:]] or [1.0,0.9,0.7,0.5,0.3]:
    Q=np.sqrt(1-a*a); w=m*a/(1+a*a); k=2*w; A=specA(a*w); lam=A+a*a*w*w-2*a*m*w
    disc=(s+.5)**2-k*k+lam
    if disc>=0: print(f"a={a} Q={Q:.3f}: real weights (delta^2={-disc:.3f}) - not oscillatory"); continue
    d=np.sqrt(-disc); pp=-(s+.5)+1j*d; pm=-(s+.5)-1j*d
    dec=lambda x,R,dR: np.linalg.solve(np.array([[x**pp,x**pm],[pp*x**(pp-1),pm*x**(pm-1)]]),[R,dR])
    def near(n,tau=1.0):
        def f(x,R,dR):
            D=x*(x+tau); U=(k*x+n*tau/2)**2-1j*s*(2*x+tau)*(k*x+n*tau/2)
            return dR,-((s+1)*(2*x+tau)/D*dR+(U/D**2+(2j*s*k-lam)/D)*R)
        q=-s-1j*n/2; x0=1e-7; a_,b_=dec(1e7,*integ(f,(x0**q,q*x0**(q-1)),x0,1e7)); return abs(b_/a_)
    def ff(r,R,dR):
        D=(r-1)**2; K=(r*r+a*a)*w-a*m
        return dR,-((s+1)*2*(r-1)/D*dR+((K*K-2j*s*(r-1)*K)/D+4j*s*w*r-lam)/D*R)
    r1=400.; rs=r1+2*np.log(r1-1)-(1+a*a)/(r1-1)
    R1=r1**(-1-2*s)*np.exp(1j*w*rs); dR1=R1*((-1-2*s)/r1+1j*w*(r1*r1+a*a)/(r1-1)**2)
    a_,b_=dec(1e-4,*integ(ff,(R1,dR1),r1,1+1e-4)); F=abs(b_/a_)
    band=[near(n) for n in (-20,-2,-.5,-.1,-.01)]
    print(f"a={a} Q={Q:.3f} delta={d:.3f}: |F|={F:.4g}  |N| band in [{min(band):.4g},{max(band):.4g}]  ratio={F/max(band):.3g}")

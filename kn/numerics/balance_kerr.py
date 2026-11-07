# Near-extremal Kerr balance test, s=-2, l=m=2.  Near zone (exact hypergeometric model) vs far zone (extremal Kerr).
import numpy as np, mpmath as mp
from scipy.integrate import solve_ivp
from sympy.physics.wigner import clebsch_gordan as CG
s,m=-2,2
# spheroidal eigenvalue at c=a*omega=m/2
ls=list(range(max(abs(s),abs(m)),m+40)); n_=len(ls); C=np.zeros((n_,n_))
for i,l in enumerate(ls):
    for j,lp in enumerate(ls):
        if abs(l-lp)<=1: C[i,j]=float(np.sqrt((2*lp+1)/(2*l+1))*CG(lp,1,l,m,0,m)*CG(lp,1,l,-s,0,-s))
c=m/2; H=np.diag([l*(l+1)-s*(s+1) for l in ls])-c*c*(C@C)+2*c*s*C
A=np.sort(np.linalg.eigvals(H).real)[0]
lam=A+c*c-2*m*c          # Teukolsky lambda with a=1, omega=m/2
d=np.sqrt(7*m*m/4-(s+.5)**2-A); pp=-(s+.5)+1j*d; pm=-(s+.5)-1j*d
print(f"A={A:.5f} lambda={lam:.5f} delta={d:.5f}")
def integ(f,y0,x0,x1):
    def F(x,y):
        z=y[0]+1j*y[1]; w=y[2]+1j*y[3]; dz,dw=f(x,z,w)
        return [dz.real,dz.imag,dw.real,dw.imag]
    r=solve_ivp(F,(x0,x1),[y0[0].real,y0[0].imag,y0[1].real,y0[1].imag],rtol=1e-11,atol=1e-14,method='DOP853')
    return r.y[0,-1]+1j*r.y[1,-1], r.y[2,-1]+1j*r.y[3,-1]
def decompose(x,R,dR):   # R = a x^pp + b x^pm
    M=np.array([[x**pp,x**pm],[pp*x**(pp-1),pm*x**(pm-1)]]); return np.linalg.solve(M,[R,dR])
# ---- near zone, tau=1
def near_ratio(n):
    tau=1.0
    def f(x,R,dR):
        D=x*(x+tau); U=(m*x+n*tau/2)**2-1j*s*(2*x+tau)*(m*x+n*tau/2)
        return dR, -((s+1)*(2*x+tau)/D*dR+(U/D**2+(2j*s*m-lam)/D)*R)
    q=-s-1j*n/2; x0=1e-7
    a,b=decompose(1e7,*integ(f,(x0**q,q*x0**(q-1)),x0,1e7))  # crude: horizon Frobenius leading term
    return b/a
# ---- far zone, extremal Kerr M=a=1, omega=m/2 ; outgoing at infinity, integrate inward
w=m/2
def ffar(r,R,dR):
    D=(r-1)**2; K=(r*r+1)*w-m
    return dR, -((s+1)*2*(r-1)/D*dR+((K*K-2j*s*(r-1)*K)/D+4j*s*w*r-lam)/D*R)
r1=400.0
rs=lambda r: r+2*np.log(r-1)-2/(r-1)        # tortoise for extremal (M=1): r*=r+2ln(r-1)-2/(r-1)
R1=r1**(-1-2*s)*np.exp(1j*w*rs(r1)); dR1=R1*((-1-2*s)/r1+1j*w*(r1*r1+1)/(r1-1)**2)
Rx,dRx=integ(ffar,(R1,dR1),r1,1+1e-4)
a,b=decompose(1e-4,Rx,dRx); F=b/a
print(f"|F| far = {abs(F):.6f}")
for n in [-20,-10,-5,-3,-2,-1,-.5,-.2,-.05]:
    print(f"n={n:6}: |N|={abs(near_ratio(n)):.6f}")

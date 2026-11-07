# At the real/oscillatory transition (delta=0) the throat basis is x^q, x^q log x (q=-(s+1/2)).
# Scaling by tau shifts 1/rho by real L=log tau, so a real mode needs Im(1/rho_N(n)) = Im(1/rho_F) exactly.
import numpy as np, sys
from scipy.optimize import brentq
exec(open('bkn_s.py').read().split('for a in')[0])
def disc_of(a):
    w=m*a/(1+a*a); k=2*w; A=specA(a*w); lam=A+a*a*w*w-2*a*m*w; return (s+.5)**2-k*k+lam
astar=brentq(disc_of,0.5,0.7); a=astar
w=m*a/(1+a*a); k=2*w; A=specA(a*w); lam=A+a*a*w*w-2*a*m*w; q=-(s+.5)
dec=lambda x,R,dR: np.linalg.solve(np.array([[x**q,x**q*np.log(x)],[q*x**(q-1),x**(q-1)*(q*np.log(x)+1)]]),[R,dR])
def invrho_near(n,tau=1.0):
    def f(x,R,dR):
        D=x*(x+tau); U=(k*x+n*tau/2)**2-1j*s*(2*x+tau)*(k*x+n*tau/2)
        return dR,-((s+1)*(2*x+tau)/D*dR+(U/D**2+(2j*s*k-lam)/D)*R)
    p=-s-1j*n/2; x0=1e-7; al,be=dec(1e7,*integ(f,(x0**p,p*x0**(p-1)),x0,1e7)); return al/be
def ff(r,R,dR):
    D=(r-1)**2; K=(r*r+a*a)*w-a*m
    return dR,-((s+1)*2*(r-1)/D*dR+((K*K-2j*s*(r-1)*K)/D+4j*s*w*r-lam)/D*R)
r1=400.; rs=r1+2*np.log(r1-1)-(1+a*a)/(r1-1)
R1=r1**(-1-2*s)*np.exp(1j*w*rs); dR1=R1*((-1-2*s)/r1+1j*w*(r1*r1+a*a)/(r1-1)**2)
al,be=dec(1e-4,*integ(ff,(R1,dR1),r1,1+1e-4)); iF=al/be
print(f"m={m} transition a*={a:.5f} Q*={np.sqrt(1-a*a):.5f}; Im(1/rho_F)={iF.imag:.5f}")
for n in [-50,-10,-2,-.5,-.1,-.01,-.001]:
    print(f"  n={n:7}: Im(1/rho_N)={invrho_near(n).imag:.5f}")

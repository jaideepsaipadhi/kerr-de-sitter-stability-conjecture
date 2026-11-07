# B1 checks: (1) R''(r_c)>0 identity at spherical orbits; (2) numerically: R has <=2 roots in (r+,inf) when Theta>=0 allowed.
import sympy as sp, numpy as np
r,M,a,Q,lam,eta=sp.symbols('r M a Q lambda eta',real=True)
D=r**2-2*M*r+a**2+Q**2; P=r**2+a**2-a*lam; R=P**2-D*eta
# spherical: P = 2 r D/(r-M), eta = P^2/D
Ps=2*r*D/(r-M); Rpp=sp.diff(R,r,2).subs(eta,Ps**2/D).subs(lam,(r**2+a**2-Ps)/a)
target=8*r/(r-M)**2*((r-M)**3+M*(M**2-a**2-Q**2))
print('R\'\' identity:',sp.simplify(Rpp-target))
# check R'=0 consistency
Rp=sp.diff(R,r).subs(eta,Ps**2/D).subs(lam,(r**2+a**2-Ps)/a); print('R\' =0:',sp.simplify(Rp))
rng=np.random.default_rng(0); bad=0; N=200000
for _ in range(N):
    th=rng.uniform(0,np.pi/2); rho=rng.uniform(0,0.999); A=rho*np.cos(th); q=rho*np.sin(th)
    rp=1+np.sqrt(1-A*A-q*q); L=rng.uniform(-30,30)
    # Theta>=0 for some theta: K >= min_{s in (0,1]} (L - A s)^2/s ; use eta=K
    s=np.linspace(1e-4,1,400); kmin=np.min((L-A*s)**2/s)
    K=kmin+rng.exponential(5)*rng.choice([0.01,1,30])
    c=np.polynomial.polynomial.Polynomial
    Pp=c([A*A-A*L,0,1]); Dp=c([A*A+q*q,-2,1]); Rp_=Pp**2-Dp*K
    rts=Rp_.roots(); ext=[z.real for z in rts if abs(z.imag)<1e-9 and z.real>rp+1e-9]
    if len(ext)>2: bad+=1; print('>2 roots',A,q,L,K,sorted(ext)) if bad<5 else None
print('violations',bad,'of',N)

# Near-horizon coupled weights for extremal KN: psi_-2 = x^h f2(th), psi_-1 = x^(h-1) f1(th).
# Galerkin in Jacobi-type spin-weighted bases; operator is quadratic in h -> companion linearization.
import sympy as sp, mpmath as mp, numpy as np, sys
from scipy.special import eval_jacobi
from numpy.polynomial.legendre import leggauss
src=open('nh_coupled.py').read().split('def lead')[0].replace("f2,f1=sp.Function('f2')(th),sp.Function('f1')(th)","")
exec(src)
av=sp.Rational(sys.argv[1]); mv=int(sys.argv[2]); N=int(sys.argv[3]) if len(sys.argv)>3 else 7
sg=int(sys.argv[4]) if len(sys.argv)>4 else 1
c=sp.cos(th)
def basis(s,n):
    al_,be_=abs(mv+sg*s),abs(mv-sg*s)
    return (1-c)**sp.Rational(al_,2)*(1+c)**sp.Rational(be_,2)*sp.jacobi(n,al_,be_,c)
B2=[basis(-2,n) for n in range(N)]; B1=[basis(-1,n) for n in range(N)]
xs,ws=leggauss(40); ths=np.arccos(xs)
X=mp.mpf('1e-18'); mp.mp.dps=60
def coeff(e,base,hv):
    f=sp.lambdify((th,),(e/x**base).subs({a:av,m:mv,h:hv,x:sp.Float(X,60)}),'mpmath')
    return np.array([complex(f(mp.mpf(t))) for t in ths])
def Mat(hv):
    blocks=np.zeros((2*N,2*N),complex)
    for j in range(2*N):
        p2=x**h*B2[j] if j<N else 0*x
        p1=0*x if j<N else x**(h-1)*B1[j-N]
        F2,G2,H2,F1,G1,H1=ops(p2,p1)
        e2=sp.expand((F2+Q**2*(G2+H2)).subs({a:av,m:mv}))
        e1=sp.expand((F1+Q**2*(G1+H1)).subs({a:av,m:mv}))
        v2=coeff(e2,h,hv); v1=coeff(e1,h-1,hv)
        for i in range(N):
            t2=np.array([complex(sp.N(B2[i].subs(th,t))) for t in ths]); t1=np.array([complex(sp.N(B1[i].subs(th,t))) for t in ths])
            blocks[i,j]=np.sum(ws*t2*v2); blocks[N+i,j]=np.sum(ws*t1*v1)
    return blocks
M0,M1,M2=Mat(0),Mat(1),Mat(2)
A2=(M2-2*M1+M0)/2; A1=M1-M0-A2; A0=M0
n=2*N; Z=np.zeros((n,n)); I=np.eye(n)
Lft=np.block([[Z,I],[-A0,-A1]]); Rgt=np.block([[I,Z],[Z,A2]])
from scipy.linalg import eig
ev=eig(Lft,Rgt,right=False); ev=ev[np.isfinite(ev)&(abs(ev)<40)]
ev=sorted(ev,key=lambda z:abs(z.real-1.5)+0.01*abs(z.imag))
print(f"a={av} Q={float(sp.sqrt(1-av**2)):.3f} m={mv} N={N} sg={sg}: h closest to Re=3/2:",[f"{z:.4f}" for z in ev[:6]])

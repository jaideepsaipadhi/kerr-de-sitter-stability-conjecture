# Fast coupled near-horizon weights (extremal KN, omega=m*Omega_H). Symbolic ops ONCE with generic f's,
# coefficients of (f,f',f'') lambdified with args (x,theta,h,a,m); Galerkin in Jacobi bases; quadratic in h.
import sympy as sp, mpmath as mp, numpy as np, sys, pickle, os
from numpy.polynomial.legendre import leggauss
from scipy.linalg import eig
exec(open('nh_coupled.py').read().split('def lead')[0])
cache='fw_cache.pkl'
if os.path.exists(cache): C=pickle.load(open(cache,'rb'))
else:
    d0,d1,d2,e0,e1,e2=sp.symbols('d0 d1 d2 e0 e1 e2')
    F2,G2,H2,F1,G1,H1=ops(x**h*f2,x**(h-1)*f1)
    def strip(e,base):
        e=(e/x**base).doit()
        e=e.subs({sp.Derivative(f2,(th,2)):d2,sp.Derivative(f2,th):d1,f2:d0,sp.Derivative(f1,(th,2)):e2,sp.Derivative(f1,th):e1,f1:e0})
        vs=(d0,d1,d2,e0,e1,e2)
        return [e.subs({w:(1 if w==v else 0) for w in vs}) for v in vs]
    C={'r2':strip(F2+Q**2*(G2+H2),h),'r1':strip(F1+Q**2*(G1+H1),h-1)}
    pickle.dump(C,open(cache,'wb'))
mp.mp.dps=50; X=mp.mpf('1e-20')
fn={k:[sp.lambdify((x,th,h,a,m),c,'mpmath') for c in v] for k,v in C.items()}
def run(av,mv,N=8,nq=48):
    xs,ws=leggauss(nq); ths=np.arccos(xs)
    def basis(s):
        al_,be_=abs(mv+s),abs(mv-s); out=[]
        for n in range(N):
            b=(1-sp.cos(th))**sp.Rational(al_,2)*(1+sp.cos(th))**sp.Rational(be_,2)*sp.jacobi(n,al_,be_,sp.cos(th))
            out.append([np.array([complex(sp.lambdify(th,sp.diff(b,th,k),'mpmath')(mp.mpf(t))) for t in ths]) for k in range(3)])
        return out
    B2,B1=basis(-2),basis(-1)
    def Mat(hv):
        coef={k:[np.array([complex(f(X,mp.mpf(t),hv,mp.mpf(str(av)),mv)) for t in ths]) for f in v] for k,v in fn.items()}
        A=np.zeros((2*N,2*N),complex)
        for j in range(2*N):
            g=[B2[j][k] if j<N else 0*ths for k in range(3)]+[B1[j-N][k] if j>=N else 0*ths for k in range(3)]
            v2=sum(coef['r2'][k]*g[k] for k in range(6)); v1=sum(coef['r1'][k]*g[k] for k in range(6))
            for i in range(N):
                A[i,j]=np.sum(ws*B2[i][0]*v2); A[N+i,j]=np.sum(ws*B1[i][0]*v1)
        return A
    M0,M1,M2=Mat(0),Mat(1),Mat(2); A2=(M2-2*M1+M0)/2; A1=M1-M0-A2; A0=M0
    n=2*N; Z=np.zeros((n,n)); I=np.eye(n)
    ev=eig(np.block([[Z,I],[-A0,-A1]]),np.block([[I,Z],[Z,A2]]),right=False)
    ev=ev[np.isfinite(ev)&(abs(ev)<30)]
    osc=sorted({round(abs(z.imag),4) for z in ev if abs(z.real-1.5)<1e-3 and abs(z.imag)>1e-6},reverse=True)
    real=sorted({round(abs(z.real-1.5),4) for z in ev if abs(z.imag)<1e-6 and z.real>1.5})
    return osc,real
if __name__=='__main__':
    mv=int(sys.argv[1])
    for av in [float(v) for v in sys.argv[2:]]:
        osc,real=run(av,mv); print(f"m={mv} a={av:.3f} Q={np.sqrt(1-av*av):.3f}: delta(osc)={osc[:3]}  nu(real)={real[:3]}",flush=True)

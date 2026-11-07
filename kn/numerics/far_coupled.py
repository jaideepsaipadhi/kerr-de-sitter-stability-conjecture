# Coupled far-zone solve on extremal KN at threshold (omega = m Omega_H) + near/far balance test, oscillatory channel.
# psi_-2 = U(r) f2(th), psi_-1 = V(r) f1(th), Galerkin in spin-weighted Jacobi bases; radial ODE system integrated inward
# from outgoing data at large r; near-throat decomposition into near-horizon channels x^h [v2; x^-1 v1].
import sympy as sp, numpy as np, sys, os, pickle
from numpy.polynomial.legendre import leggauss
from scipy.integrate import solve_ivp
from scipy.linalg import eig, null_space
exec(open('nh_coupled.py').read().split('def lead')[0].replace("f2,f1=sp.Function('f2')(th),sp.Function('f1')(th)",""))
U,V=sp.Function('U')(x),sp.Function('V')(x); g2,g1=sp.Function('g2')(th),sp.Function('g1')(th)
cache='far_cache.pkl'
if os.path.exists(cache): CO=pickle.load(open(cache,'rb'))
else:
    F2,G2,H2,F1,G1,H1=ops(U*g2,V*g1)
    syms=sp.symbols('u0 u1 u2 v0 v1 v2 a0 a1 a2 b0 b1 b2')
    rep={sp.Derivative(U,(x,2)):syms[2],sp.Derivative(U,x):syms[1],sp.Derivative(V,(x,2)):syms[5],sp.Derivative(V,x):syms[4],
         sp.Derivative(g2,(th,2)):syms[8],sp.Derivative(g2,th):syms[7],sp.Derivative(g1,(th,2)):syms[11],sp.Derivative(g1,th):syms[10]}
    CO={}
    for name,e in (('r2',F2+Q**2*(G2+H2)),('r1',F1+Q**2*(G1+H1))):
        e=e.doit().subs(rep).subs({U:syms[0],V:syms[3],g2:syms[6],g1:syms[9]}); e=sp.expand(e)
        # bilinear pieces: (U^k * g2^l) and (V^k * g1^l)
        d={}
        for k in range(3):
            for l in range(3):
                d[('U',k,l)]=e.coeff(syms[k]).coeff(syms[6+l]); d[('V',k,l)]=e.coeff(syms[3+k]).coeff(syms[9+l])
        CO[name]=d
    pickle.dump(CO,open(cache,'wb'))
def run(av,mv,N=8,nq=48,x0=2e-3,x1=400.):
    fn={nm:{key:sp.lambdify((x,th),c.subs({a:av,m:mv}),'numpy') for key,c in d.items()} for nm,d in CO.items()}
    xs,ws=leggauss(nq); ths=np.arccos(xs)
    def basis(s):
        al_,be_=abs(mv+s),abs(mv-s); out=[]
        for n in range(N):
            b=(1-sp.cos(th))**sp.Rational(al_,2)*(1+sp.cos(th))**sp.Rational(be_,2)*sp.jacobi(n,al_,be_,sp.cos(th))
            out.append(np.array([sp.lambdify(th,sp.diff(b,th,k),'numpy')(ths)*np.ones(nq) for k in range(3)]))
        return np.array(out)  # N x 3 x nq
    B2,B1=basis(-2),basis(-1)
    def ev(f,X): v=f(X,ths); return v*np.ones(nq)
    def mats(X):
        A=[np.zeros((2*N,2*N),complex) for _ in range(3)]
        for (eq,T) in (('r2',B2),('r1',B1)):
            row=0 if eq=='r2' else N
            for (blk,Bs,col) in (('U',B2,0),('V',B1,N)):
                for k in range(3):
                    c=[ev(fn[eq][(blk,k,l)],X) for l in range(3)]
                    act=sum(c[l][None,:]*Bs[:,l,:] for l in range(3))      # N(j) x nq
                    A[k][row:row+N,col:col+N]=np.einsum('q,iq,jq->ij',ws,T[:,0,:],act)
        return A
    w=mv*av/(1+av*av); Q2=1-av*av
    def rhs(X,y):
        A0,A1,A2=mats(X); n2=2*N
        Y=y[:n2*2*n2].view(complex).reshape(n2,2*n2) if False else None
        return None
    # integrate the 2N outgoing solutions as a matrix ODE
    n2=2*N
    def F(X,y):
        Z=y.view(complex).reshape(2,n2,n2); A0,A1,A2=mats(X)
        d2=-np.linalg.solve(A2,A1@Z[1]+A0@Z[0]); return np.stack([Z[1],d2]).ravel().view(float)
    rs=lambda X:(1+X)+2*np.log(X)-(1+av*av)/X
    p=np.array([3.0]*N+[1.0]*N)       # outgoing powers r^(-1-2s): s=-2 -> 3, s=-1 -> 1
    r1=1+x1; ph=np.exp(1j*w*rs(x1))*r1**p
    dph=ph*(p/r1+1j*w*(r1*r1+av*av)/x1**2)
    Z0=np.stack([np.diag(ph),np.diag(dph)])
    sol=solve_ivp(F,(x1,x0),Z0.ravel().view(float),rtol=1e-8,atol=1e-12,method='DOP853')
    Zx=sol.y[:,-1].view(complex).reshape(2,n2,n2)
    # near-horizon channels at extremality from the quadratic eigenproblem (reuse fast_weights machinery)
    import importlib.util; spec=importlib.util.spec_from_file_location('fw','fast_weights.py'); fw=importlib.util.module_from_spec(spec)
    sys.argv=['fw']; spec.loader.exec_module(fw)
    return Zx, w

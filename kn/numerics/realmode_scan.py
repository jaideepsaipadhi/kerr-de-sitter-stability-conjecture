# Direct real-mode test for the FULL coupled DGS system on near-extremal Kerr-Newman.
# For real omega in the superradiant band: 2N horizon-regular solutions (outward) vs 2N outgoing (inward);
# smallest singular value of the column-normalised matching matrix -> 0 iff a real mode exists.
import sympy as sp, numpy as np, sys, os, pickle, time
from numpy.polynomial.legendre import leggauss
from scipy.integrate import solve_ivp
r,th,a,Q,m,w=sp.symbols('r theta a Q m omega'); M=1
Dl=r**2-2*M*r+a**2+Q**2
Kr=a*m-(r**2+a**2)*w; Kt=m/sp.sin(th)-a*w*sp.sin(th)
rb=r+sp.I*a*sp.cos(th); rbs=r-sp.I*a*sp.cos(th)
al=lambda sg: 1/(3*(rb**2*M-rb*Q**2)+sg*Q**2*rbs)
D =lambda j,u: sp.diff(u,r)+sp.I*Kr/Dl*u+2*j*(r-M)/Dl*u
Dd=lambda j,u: sp.diff(u,r)-sp.I*Kr/Dl*u+2*j*(r-M)/Dl*u
L =lambda j,u: sp.diff(u,th)+Kt*u+j*sp.cot(th)*u
Ld=lambda j,u: sp.diff(u,th)-Kt*u+j*sp.cot(th)*u
ia=sp.I*a*sp.sin(th)
U,V=sp.Function('U')(r),sp.Function('V')(r); g2,g1=sp.Function('g2')(th),sp.Function('g1')(th)
cache='rm_cache.pkl'
if os.path.exists(cache): CO=pickle.load(open(cache,'rb'))
else:
    p2,p1=U*g2,V*g1
    F2=Dl*Dd(-1,D(0,p2))+L(-1,Ld(2,p2))-6*sp.I*w*rb*p2
    G2=Dl*Dd(-1,al(-1)*rbs*D(0,p2))-3*Dl*Dd(-1,al(-1)*p2)-L(-1,al(1)*rbs*Ld(2,p2))+3*L(-1,al(1)*ia*p2)
    H2=-Dl*Dd(-1,al(-1)*rbs*L(-1,p1))-3*Dl*Dd(-1,al(-1)*ia*p1)-L(-1,al(1)*rbs*Dl*Dd(-1,p1))-3*L(-1,al(1)*Dl*p1)
    F1=Dl*D(1,Dd(-1,p1))+Ld(2,L(-1,p1))-6*sp.I*w*rb*p1
    G1=-D(0,al(1)*rbs*Dl*Dd(-1,p1))-3*D(0,al(1)*Dl*p1)+Ld(2,al(-1)*rbs*L(-1,p1))+3*Ld(2,al(-1)*ia*p1)
    H1=-D(0,al(1)*rbs*Ld(2,p2))+3*D(0,al(1)*ia*p2)-Ld(2,al(-1)*rbs*D(0,p2))+3*Ld(2,al(-1)*p2)
    s_=sp.symbols('u0 u1 u2 v0 v1 v2 a0 a1 a2 b0 b1 b2')
    rep={sp.Derivative(U,(r,2)):s_[2],sp.Derivative(U,r):s_[1],sp.Derivative(V,(r,2)):s_[5],sp.Derivative(V,r):s_[4],
         sp.Derivative(g2,(th,2)):s_[8],sp.Derivative(g2,th):s_[7],sp.Derivative(g1,(th,2)):s_[11],sp.Derivative(g1,th):s_[10]}
    CO={}
    for name,e in (('r2',F2+Q**2*(G2+H2)),('r1',F1+Q**2*(G1+H1))):
        e=sp.expand(e.doit().subs(rep).subs({U:s_[0],V:s_[3],g2:s_[6],g1:s_[9]}))
        z={v:0 for v in s_}
        for k in range(3):
            for l in range(3):
                CO[(name,'U',k,l)]=e.subs({**z,s_[k]:1,s_[6+l]:1})
                CO[(name,'V',k,l)]=e.subs({**z,s_[3+k]:1,s_[9+l]:1})
    pickle.dump(CO,open(cache,'wb'))
def setup(av,Qv,mv,N=8,nq=40):
    fn={k:sp.lambdify((r,th,w),c.subs({a:av,Q:Qv,m:mv}),'numpy') for k,c in CO.items()}
    xs,ws=leggauss(nq); ths=np.arccos(xs)
    def basis(s):
        A_,B_=abs(mv+s),abs(mv-s); out=[]
        for n in range(N):
            b=(1-sp.cos(th))**sp.Rational(A_,2)*(1+sp.cos(th))**sp.Rational(B_,2)*sp.jacobi(n,A_,B_,sp.cos(th))
            out.append(np.array([sp.lambdify(th,sp.diff(b,th,k),'numpy')(ths)*np.ones(nq) for k in range(3)]))
        return np.array(out)
    B2,B1=basis(-2),basis(-1)
    def mats(R,W):
        A=[np.zeros((2*N,2*N),complex) for _ in range(3)]
        for eq,T,row in (('r2',B2,0),('r1',B1,N)):
            for blk,Bs,col in (('U',B2,0),('V',B1,N)):
                for k in range(3):
                    act=sum((fn[(eq,blk,k,l)](R,ths,W)*np.ones(nq))[None,:]*Bs[:,l,:] for l in range(3))
                    A[k][row:row+N,col:col+N]=np.einsum('q,iq,jq->ij',ws,T[:,0,:],act)
        return A
    return mats
def smin(av,Qv,mv,W,mats,N=8,x0=1e-4,rmid=None,rbig=300.):
    rp=1+np.sqrt(1-av**2-Qv**2); rm_=1-np.sqrt(1-av**2-Qv**2); OH=av/(rp**2+av**2)
    sig=(W-mv*OH)*(rp**2+av**2)/(rp-rm_); n2=2*N
    def F(R,y):
        Z=y.view(complex).reshape(2,n2,n2); A0,A1,A2=mats(R,W)
        return np.stack([Z[1],-np.linalg.solve(A2,A1@Z[1]+A0@Z[0])]).ravel().view(float)
    q=np.array([2-1j*sig]*N+[1-1j*sig]*N); X=x0
    Zh=np.stack([np.diag(X**q),np.diag(q*X**(q-1))])
    rmid=rmid or rp+1.0
    a1=solve_ivp(F,(rp+X,rmid),Zh.ravel().view(float),rtol=1e-8,atol=1e-30,method='DOP853').y[:,-1].view(complex).reshape(2,n2,n2)
    rs=lambda R:R+(rp**2+av**2)/(rp-rm_)*np.log(R-rp)-(rm_**2+av**2)/(rp-rm_)*np.log(R-rm_)
    p=np.array([3.]*N+[1.]*N); ph=np.exp(1j*W*rs(rbig))*rbig**p
    dph=ph*(p/rbig+1j*W*(rbig**2+av**2)/((rbig-rp)*(rbig-rm_)))
    a2=solve_ivp(F,(rbig,rmid),np.stack([np.diag(ph),np.diag(dph)]).ravel().view(float),rtol=1e-8,atol=1e-30,method='DOP853').y[:,-1].view(complex).reshape(2,n2,n2)
    Mm=np.block([[a1[0],a2[0]],[a1[1],a2[1]]]); Mm=Mm/np.linalg.norm(Mm,axis=0)
    return np.linalg.svd(Mm,compute_uv=False)[-1]
if __name__=='__main__':
    av,eps,mv=float(sys.argv[1]),float(sys.argv[2]),int(sys.argv[3]); Qv=np.sqrt(1-eps-av*av)
    mats=setup(av,Qv,mv); rp=1+np.sqrt(eps); OH=av/(rp**2+av**2)
    for fr in [0.05,0.2,0.4,0.6,0.8,0.9,0.95,0.98,0.995]:
        W=fr*mv*OH; t=time.time()
        print(f"a={av} Q={Qv:.4f} eps={eps} m={mv} omega/(m OmH)={fr}: smin={smin(av,Qv,mv,W,mats):.3e} ({time.time()-t:.0f}s)",flush=True)

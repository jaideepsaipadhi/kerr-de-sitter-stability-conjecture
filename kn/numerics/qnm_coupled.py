# Spectral QNM solver for the FULL coupled DGS system on Kerr-Newman (Chebyshev in u=1-r+/r  x  Galerkin in theta).
# psi_s = P_s(r) h_s(u,theta), P_s=(r-r+)^(-s-i sig)(r-r-)^(-1-s+i sig+2iMw)e^{iwr}; rows divided by P_-2 (eq r2) / P_-1 (eq r1),
# then row-scaled (omega-independent). Operator polynomial in w of degree <= 4 -> companion linearisation.
import sympy as sp, numpy as np, pickle, sys, time
from numpy.polynomial.legendre import leggauss
from numpy.polynomial import chebyshev as Ch
from scipy.linalg import eig
r,th,a,Q,m,w=sp.symbols('r theta a Q m omega')
CO=pickle.load(open('rm_cache.pkl','rb'))
def build(av,Qv,mv,Nr=24,Nt=8,nq=32,M=1.0):
    rp=M+np.sqrt(M*M-av*av-Qv*Qv); rm=M-np.sqrt(M*M-av*av-Qv*Qv); OH=av/(rp*rp+av*av)
    fn={k:sp.lambdify((r,th,w),c.subs({a:av,Q:Qv,m:mv}),'numpy') for k,c in CO.items()}
    k=np.arange(Nr); x=np.cos(np.pi*(2*k+1)/(2*Nr)); uu=(x+1)/2; rr=rp/(1-uu)
    T=[np.array([Ch.chebval(x,Ch.chebder(np.eye(Nr)[n],d))*2**d for n in range(Nr)]).T for d in range(3)]  # node x n
    ur=(1-uu)**2/rp; urr=-2*(1-uu)**3/rp**2
    xq,wq=leggauss(nq); tq=np.arccos(xq)
    def basis(s):
        A_,B_=abs(mv+s),abs(mv-s); out=[]
        for n in range(Nt):
            b=(1-sp.cos(th))**sp.Rational(A_,2)*(1+sp.cos(th))**sp.Rational(B_,2)*sp.jacobi(n,A_,B_,sp.cos(th))
            out.append(np.array([sp.lambdify(th,sp.diff(b,th,d),'numpy')(tq)*np.ones(nq) for d in range(3)]))
        return np.array(out)
    Bs={'U':basis(-2),'V':basis(-1)}; sp_={'U':-2,'V':-1}
    def lnP(s,W):
        sig=(W-mv*OH)*(rp*rp+av*av)/(rp-rm); A=-s-1j*sig; B=-1-s+1j*sig+2j*M*W
        p1=A/(rr-rp)+B/(rr-rm)+1j*W; dp1=-A/(rr-rp)**2-B/(rr-rm)**2
        lP=A*np.log(rr-rp)+B*np.log(rr-rm)+1j*W*rr
        return p1,dp1,lP
    def Mat(W):
        n=Nr*Nt; A=np.zeros((2*n,2*n),complex)
        for eq,rowoff,srow in (('r2',0,-2),('r1',n,-1)):
            _,_,lProw=lnP(srow,W)
            for blk,coloff in (('U',0),('V',n)):
                p1,dp1,lP=lnP(sp_[blk],W); ratio=np.exp(lP-lProw)
                R=[T[0], T[1]*ur[:,None]+p1[:,None]*T[0],
                   T[2]*(ur**2)[:,None]+T[1]*(urr+2*p1*ur)[:,None]+(dp1+p1**2)[:,None]*T[0]]
                for kk in range(3):
                    for l in range(3):
                        C=fn[(eq,blk,kk,l)](rr[:,None],tq[None,:],W)*np.ones((Nr,nq))     # node x quad
                        G=np.einsum('q,iq,jq,aq->aij',wq,Bs['U' if eq=='r2' else 'V'][:,0,:],Bs[blk][:,l,:],C)  # node x test x basis
                        A[rowoff:rowoff+n,coloff:coloff+n]+=np.einsum('a,aij,an->aijn',ratio,G,R[kk]).reshape(Nr,Nt,Nt*Nr).reshape(n,n) if False else \
                            np.einsum('a,aij,an->aijn',ratio,G,R[kk]).transpose(0,1,3,2).reshape(n,Nr*Nt)[:, :]*0+np.einsum('a,aij,an->ainj',ratio,G,R[kk]).reshape(n,n)
        return A
    return Mat,(Nr,Nt)
def qnms(av,Qv,mv,Nr=24,Nt=8,deg=4,wmax=2.0):
    Mat,_=build(av,Qv,mv,Nr,Nt)
    S=[0,1,-1,2,-2,0.5][:deg+1+1]; Ms=[Mat(s) for s in S]
    scale=1/np.max(np.abs(Ms[1]),axis=1)[:,None]; Ms=[scale*X for X in Ms]
    V=np.vander(np.array(S[:deg+1],complex),deg+1,increasing=True); coef=np.linalg.solve(V,np.array(Ms[:deg+1]).reshape(deg+1,-1))
    As=[c.reshape(Ms[0].shape) for c in coef]
    pred=sum(As[j]*S[deg+1]**j for j in range(deg+1)); print(" deg chk",np.abs(pred-Ms[deg+1]).max())
    n=As[0].shape[0]; I=np.eye(n); Z=np.zeros((n,n))
    L=np.zeros((deg*n,deg*n),complex); Rm=np.eye(deg*n,dtype=complex)
    for j in range(deg-1): L[j*n:(j+1)*n,(j+1)*n:(j+2)*n]=I
    for j in range(deg): L[(deg-1)*n:,j*n:(j+1)*n]=-As[j]
    Rm[(deg-1)*n:,(deg-1)*n:]=As[deg]
    ev=eig(L,Rm,right=False); ev=ev[np.isfinite(ev)&(abs(ev)<wmax)&(ev.real>0.05)]
    return sorted(ev,key=lambda z:-z.imag)
if __name__=='__main__':
    av,Qv,mv=float(sys.argv[1]),float(sys.argv[2]),int(sys.argv[3]); Nr=int(sys.argv[4]); Nt=int(sys.argv[5])
    t=time.time(); e=qnms(av,Qv,mv,Nr,Nt); print(f"a={av} Q={Qv} m={mv} Nr={Nr} Nt={Nt}:",[f"{z:.5f}" for z in e[:6]],f"{time.time()-t:.0f}s")

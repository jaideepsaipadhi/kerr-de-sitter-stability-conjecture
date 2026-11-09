import numpy as np, time, sys
from iv import I
def Dbox(t,ze,y,lr,be):
    # variables as intervals; A = r*B, r = 10**lr ; B = be*(5.5t-1.25)
    r=I(10**lr.lo,10**lr.hi)
    B=be*(5.5*t-1.25); A=r*B
    x0=ze*(t-1); z0=1+x0; x=y*x0; z=1+x; s=(1-y)*x0
    zi=z.inv()
    mu=z*x*(t-z); mup=-3*z.sq()+2*(1+t)*z-t
    L=A*zi+B-3*z; Lp=-(A*zi.sq())-3; Lpp=2*A*zi*zi.sq()
    Li=L.inv(); g=-(Lp*Li); gp=g.sq()-Lpp*Li
    mu0=z0*x0*(t-z0); L0=A/z0+B-3*z0
    kap=L0/mu0; w=8*kap*s.sq(); ow=(1+w).inv()
    psi=1-ow; psip=-16*kap*s*ow.sq()
    q=0.45*(1-y)/(x*z); qp=0.45*(-(1/x0)/(x*z)-(1-y)*(2*z-1)/(x*z).sq())
    G=1.7*g; h=G*psi+q; hp=1.7*gp*psi+G*psip+qp
    D=0.5*(mup*h+mu*(hp+h.sq()))-s*(L*h+Lp)
    # R3 membership: p_param>1 and 2x0^2(A/z0+3)>=t-1 ; box certainly outside if either certainly fails
    pp=2*x0*(A+B-3)/(t-1); c2=2*x0.sq()*(A/z0+3)-(t-1)
    outside=(pp.hi<=1)|(c2.hi<0)
    return D,outside
def run(lo,hi,tmax=50,cap=2_000_000):
    t0=time.time(); boxes=np.array([lo+hi],float).reshape(1,10); done=0; nb=0
    while len(boxes) and time.time()-t0<tmax:
        L_=boxes[:,:5]; H_=boxes[:,5:]
        iv=[I(L_[:,k],H_[:,k]) for k in range(5)]
        D,out=Dbox(*iv)
        ok=(D.hi<0)|out
        done+=ok.sum(); nb+=len(boxes)
        rest=boxes[~ok]
        if len(rest)==0: boxes=rest; break
        # split widest (relative) dim
        W=rest[:,5:]-rest[:,:5]; k=np.argmax(W/ (np.array(hi)-np.array(lo)),axis=1)
        mid=(rest[np.arange(len(rest)),k]+rest[np.arange(len(rest)),5+k])/2
        a=rest.copy(); b=rest.copy()
        a[np.arange(len(rest)),5+k]=mid; b[np.arange(len(rest)),k]=mid
        boxes=np.vstack([a,b])
        if len(boxes)>cap: break
    return len(boxes)==0, done, len(boxes), round(time.time()-t0,1)
if __name__=='__main__':
    lo=[2,0.2,0.02,np.log10(0.5),1]; hi=[3,0.8,0.98,np.log10(50),3]
    print(run(lo,hi))

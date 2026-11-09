# B3 certificate with mean-value form: F(box) ⊂ F(center) + sum_k dF/dx_k(box) * (box_k - c_k)
import numpy as np, time, sys
from iv import I
from di import DI, dhull, hull
def B3(t,ze,y,be,lam, sel):
    B=be*(5.5*t-1.25)
    x0=ze*(t-1); z0=1+x0; x=y*x0; z=1+x; s=(1-y)*x0
    L=(B-3*z)/lam; A=z*(L-B+3*z)
    mu=z*x*(t-z); mup=-3*z.sq()+2*(1+t)*z-t
    zi=z.inv(); N=A*zi.sq()+3; g=N/L
    v=1-y; q=0.45*v/(x*z); qp=0.45*(-(1/x0)/(x*z)-v*(2*z-1)/(x*z).sq())
    hup=1.7*g+q; hlo=(17/11)*g+q
    if sel==1: t1=0.5*mup*hup
    elif sel==-1: t1=0.5*mup*hlo
    else:
        a_,b_=0.5*mup*hup,0.5*mup*hlo
        t1=dhull(a_,b_) if isinstance(a_,DI) else hull(a_,b_)
    F=t1+0.85*mu*g.sq()+0.5*mu*hup.sq()+0.5*mu*qp-0.5455*s*L*g-s*L*q
    mu0=z0*x0*(t-z0); L0=A/z0+B-3*z0; w=8*L0*s.sq()/mu0
    pp=2*x0*(A+B-3)/(t-1); c2=2*x0.sq()*(A/z0+3)-(t-1)
    return F,w,pp,c2,mup
def val(x): return x.v if isinstance(x,DI) else x
def evalbox(boxes):
    n=len(boxes); lo=boxes[:,:5]; hi=boxes[:,5:]; c=(lo+hi)/2
    X=[DI(I(lo[:,k],hi[:,k]),[I.c(np.full(n,1.0 if j==k else 0.0)) for j in range(5)]) for k in range(5)]
    Xc=[I.c(c[:,k]) for k in range(5)]
    # sign of mu' on the box (interval) decides the branch
    _,_,_,_,mupb=B3(*[x.v for x in X],sel=1)
    sel=np.where(mupb.lo>=0,1,np.where(mupb.hi<=0,-1,0))
    out_F_hi=np.full(n,np.inf); outside=np.zeros(n,bool)
    for sv in (1,-1,0):
        m=sel==sv
        if not m.any(): continue
        Xs=[DI(I(x.v.lo[m],x.v.hi[m]),[I(gg.lo[m],gg.hi[m]) for gg in x.g]) for x in X]
        F,w,pp,c2,_=B3(*Xs,sel=sv)
        Fc,_,_,_,_=B3(*[I(xc.lo[m],xc.hi[m]) for xc in Xc],sel=sv)
        mv=Fc.hi.copy()
        if sv==0: mv=np.full(m.sum(),np.inf)
        for k in range(5 if sv!=0 else 0):
            dk=I(lo[m,k]-c[m,k],hi[m,k]-c[m,k]); term=F.g[k]*dk; mv=np.nextafter(mv+term.hi,np.inf)
        out_F_hi[m]=np.minimum(F.v.hi,mv)
        outside[m]=(w.v.hi<10)|(pp.v.hi<=1)|(c2.v.hi<0)
    return out_F_hi,outside
def run(lo,hi,tmax=50,cap=2_000_000):
    t0=time.time(); lo=np.array(lo,float); hi=np.array(hi,float)
    boxes=np.concatenate([lo,hi])[None,:]; done=0
    while len(boxes) and time.time()-t0<tmax:
        Fh,out=evalbox(boxes); ok=(Fh<0)|out
        done+=int(ok.sum()); rest=boxes[~ok]
        if len(rest)==0: boxes=rest; break
        W=(rest[:,5:]-rest[:,:5])/(hi-lo); k=np.argmax(W,axis=1); ix=np.arange(len(rest))
        mid=(rest[ix,k]+rest[ix,5+k])/2; a=rest.copy(); b=rest.copy(); a[ix,5+k]=mid; b[ix,k]=mid
        boxes=np.vstack([a,b])
        if len(boxes)>cap: break
    return dict(certified=len(boxes)==0,closed=done,open=len(boxes),sec=round(time.time()-t0,1))
if __name__=='__main__':
    a=[float(v) for v in sys.argv[1:11]]
    print(run(a[:5],a[5:]),flush=True)

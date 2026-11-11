# Certificate attempt for the bulk inequality B3 (item 107) by interval branch-and-bound.
# Variables: t, ze (zf), y, be = B/floor in [1, Bmax], lam = (B-3z)/L in (0,1]  (L in [B-3z, inf)).
import numpy as np, time, sys
from iv import I
def B3box(t,ze,y,be,lam):
    B=be*(5.5*t-1.25)
    x0=ze*(t-1); z0=1+x0; x=y*x0; z=1+x; s=(1-y)*x0
    L=(B-3*z)/lam; A=z*(L-B+3*z)
    mu=z*x*(t-z); mup=-3*z.sq()+2*(1+t)*z-t
    zi=z.inv(); N=A*zi.sq()+3; g=N/L
    v=1-y; q=0.45*v/(x*z); qp=0.45*(-(1/x0)/(x*z)-v*(2*z-1)/(x*z).sq())
    hup=1.7*g+q; hlo=(17/11)*g+q
    # 1/2 mu' h_pm: if mu'>0 use hup, if <0 use hlo; for straddling boxes take hull of both products
    t1a=0.5*mup*hup; t1b=0.5*mup*hlo
    lo=np.where(mup.lo>=0,t1a.lo,np.where(mup.hi<=0,t1b.lo,np.minimum(t1a.lo,t1b.lo)))
    hi=np.where(mup.lo>=0,t1a.hi,np.where(mup.hi<=0,t1b.hi,np.maximum(t1a.hi,t1b.hi)))
    t1=I(lo,hi)
    lhs=t1+0.85*mu*g.sq()+0.5*mu*hup.sq()+0.5*mu*qp
    rhs=0.5455*s*L*g+s*L*q
    F=lhs-rhs
    mu0=z0*x0*(t-z0); L0=A/z0+B-3*z0; w=8*L0*s.sq()/mu0
    pp=2*x0*(A+B-3)/(t-1); c2=2*x0.sq()*(A/z0+3)-(t-1)
    out=(w.hi<10)|(pp.hi<=1)|(c2.hi<0)
    return F,out
def run(lo,hi,tmax=50,cap=3_000_000):
    t0=time.time(); lo=np.array(lo,float); hi=np.array(hi,float)
    boxes=np.concatenate([lo,hi])[None,:]; done=0
    while len(boxes) and time.time()-t0<tmax:
        iv=[I(boxes[:,k],boxes[:,5+k]) for k in range(5)]
        F,out=B3box(*iv); ok=(F.hi<0)|out
        done+=int(ok.sum()); rest=boxes[~ok]
        if len(rest)==0: boxes=rest; break
        W=(rest[:,5:]-rest[:,:5])/(hi-lo); k=np.argmax(W,axis=1); ix=np.arange(len(rest))
        mid=(rest[ix,k]+rest[ix,5+k])/2; a=rest.copy(); b=rest.copy(); a[ix,5+k]=mid; b[ix,k]=mid
        boxes=np.vstack([a,b])
        if len(boxes)>cap: break
    return dict(certified=len(boxes)==0,closed=done,open=len(boxes),sec=round(time.time()-t0,1))
if __name__=='__main__':
    tlo,thi=float(sys.argv[1]),float(sys.argv[2])
    print(run([tlo,0.02,0.02,1,1e-3],[thi,0.98,0.98,3,1.0]))

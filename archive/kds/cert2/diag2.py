import numpy as np, time
from iv import I
import b3mv
lo=np.array([2,0.2,0.1,1,0.1]); hi=np.array([3,0.8,0.9,3,1.0])
boxes=np.concatenate([lo,hi])[None,:]
for it in range(18):
    Fh,out=b3mv.evalbox(boxes); ok=(Fh<0)|out; rest=boxes[~ok]
    W=(rest[:,5:]-rest[:,:5])/(hi-lo); k=np.argmax(W,axis=1); ix=np.arange(len(rest))
    mid=(rest[ix,k]+rest[ix,5+k])/2; a=rest.copy(); b=rest.copy(); a[ix,5+k]=mid; b[ix,k]=mid
    boxes=np.vstack([a,b])
c=(boxes[:,:5]+boxes[:,5:])/2
F,w,pp,c2,mup=b3mv.B3(*[I.c(c[:,k]) for k in range(5)],sel=1)
F2,_,_,_,_=b3mv.B3(*[I.c(c[:,k]) for k in range(5)],sel=-1)
Fc=np.where(mup.lo>=0,F.lo,F2.lo)
inb=(w.lo>=10)&(pp.lo>1)&(c2.lo>=0)
print('open',len(boxes),'centres in domain',inb.sum(),'violations F>0 in domain',(inb&(Fc>0)).sum(),'centres outside',(~inb).sum())
print('centres outside: w range',np.percentile(w.lo[~inb],[0,50,100]))

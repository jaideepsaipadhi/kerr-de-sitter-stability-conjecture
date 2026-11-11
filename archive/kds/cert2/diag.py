import numpy as np, time
from iv import I
from b3cert import B3box
lo=np.array([2,0.02,0.02,1,1e-3]); hi=np.array([3,0.98,0.98,3,1.0])
boxes=np.concatenate([lo,hi])[None,:]
for it in range(22):
    iv=[I(boxes[:,k],boxes[:,5+k]) for k in range(5)]
    F,out=B3box(*iv); ok=(F.hi<0)|out; rest=boxes[~ok]
    W=(rest[:,5:]-rest[:,:5])/(hi-lo); k=np.argmax(W,axis=1); ix=np.arange(len(rest))
    mid=(rest[ix,k]+rest[ix,5+k])/2; a=rest.copy(); b=rest.copy(); a[ix,5+k]=mid; b[ix,k]=mid
    boxes=np.vstack([a,b])
c=(boxes[:,:5]+boxes[:,5:])/2
names=['t','zf','y','be','lam']
for j in range(5):
    h,e=np.histogram(c[:,j],bins=8,range=(lo[j],hi[j])); print(names[j],h)

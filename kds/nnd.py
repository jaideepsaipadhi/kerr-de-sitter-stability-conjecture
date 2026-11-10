# n-D nonneg-Bernstein certificate of p>0 on the open unit box (zero coeffs allowed, face rule), exact.
import itertools, numpy as np, sympy as sp, time
from fractions import Fraction as Fr
from math import comb, lcm
def to_bern(P,V):
    P=sp.Poly(P,*V); deg=[P.degree(v) for v in V]
    A=np.empty([d+1 for d in deg],dtype=object); A.fill(Fr(0))
    for m,c in zip(P.monoms(),P.coeffs()): A[m]=Fr(int(sp.numer(c)),int(sp.denom(c)))
    for k,n in enumerate(deg):
        A=np.moveaxis(A,k,0); B=np.empty_like(A)
        for i in range(n+1):
            acc=A[0]*0
            for j in range(i+1): acc=acc+A[j]*Fr(comb(i,j),comb(n,j))
            B[i]=acc
        A=np.moveaxis(B,0,k)
    L=1
    for v in A.flat: L=lcm(L,v.denominator)
    return np.vectorize(lambda v:int(v*L),otypes=[object])(A)
def split(A,k):
    A=np.moveaxis(A,k,0); n=A.shape[0]-1; lv=[[A[i] for i in range(n+1)]]
    for r in range(1,n+1):
        p=lv[-1]; lv.append([p[i]+p[i+1] for i in range(len(p)-1)])
    Lc=np.stack([lv[r][0]*(2**(n-r)) for r in range(n+1)]); Rc=np.stack([lv[n-r][r]*(2**r) for r in range(n+1)])
    return np.moveaxis(Lc,0,k),np.moveaxis(Rc,0,k)
def ok_box(A,bx):
    if min(A.flat)<0: return False
    d=A.ndim
    for sel in itertools.product((0,1,2),repeat=d):
        if any((s==0 and bx[k][0]==0) or (s==1 and bx[k][1]==1) for k,s in enumerate(sel)): continue
        sl=tuple(0 if s==0 else (A.shape[k]-1 if s==1 else slice(None)) for k,s in enumerate(sel))
        sub=A[sl]; m=max(sub.flat) if isinstance(sub,np.ndarray) else sub
        if not m>0: return False
    return True
def cert(P,V,tmax=50,floor=Fr(1,2**14)):
    d=len(V); st=[(to_bern(P,V),((Fr(0),Fr(1)),)*d)]; closed=0; t0=time.time()
    while st:
        if time.time()-t0>tmax: return dict(ok=False,timeout=True,closed=closed,open=len(st))
        A,bx=st.pop()
        if ok_box(A,bx): closed+=1; continue
        if min(b[1]-b[0] for b in bx)<floor: return dict(ok=False,stuck=[(float(a),float(b)) for a,b in bx])
        k=max(range(d),key=lambda i:bx[i][1]-bx[i][0]); a,c=bx[k]; m=(a+c)/2
        L,R=split(A,k); bl=list(bx); bl[k]=(a,m); br=list(bx); br[k]=(m,c)
        st+=[(L,tuple(bl)),(R,tuple(br))]
    return dict(ok=True,closed=closed)

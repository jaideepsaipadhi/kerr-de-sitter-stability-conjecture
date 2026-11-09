import numpy as np, sympy as sp, time
from fractions import Fraction
from math import comb, lcm
def to_bern(P,V):
    deg=[P.degree(v) for v in V]
    A=np.empty([d+1 for d in deg],dtype=object); A.fill(Fraction(0))
    for mon,c in zip(P.monoms(),P.coeffs()): A[mon]=Fraction(int(sp.numer(c)),int(sp.denom(c)))
    for k,n in enumerate(deg):
        A=np.moveaxis(A,k,0); Bc=np.empty_like(A)
        for i in range(n+1):
            acc=A[0]*0
            for j in range(i+1): acc=acc+A[j]*Fraction(comb(i,j),comb(n,j))
            Bc[i]=acc
        A=np.moveaxis(Bc,0,k)
    Lm=1
    for v in A.flat: Lm=lcm(Lm,v.denominator)
    return np.vectorize(lambda v:int(v*Lm),otypes=[object])(A)
def split(A,k):
    A=np.moveaxis(A,k,0); n=A.shape[0]-1
    lv=[[A[i] for i in range(n+1)]]
    for r in range(1,n+1):
        p=lv[-1]; lv.append([p[i]+p[i+1] for i in range(len(p)-1)])
    Lc=np.stack([lv[r][0]*(2**(n-r)) for r in range(n+1)])
    Rc=np.stack([lv[n-r][r]*(2**r) for r in range(n+1)])
    return np.moveaxis(Lc,0,k), np.moveaxis(Rc,0,k)
def mn(A): return min(A.flat)
def mx(A): return max(A.flat)
def bb(BB,rule,ndim,tmax=48,maxdep=80,box=None):
    box=box or tuple((0.0,1.0) for _ in range(ndim))
    start=time.time(); stack=[(BB,0,box)]; closed=0; deep=[]
    while stack and time.time()-start<tmax:
        b,dep,bx=stack.pop()
        if rule(b): closed+=1; continue
        if dep>=maxdep: deep.append(bx); continue
        k=max(range(ndim),key=lambda i:bx[i][1]-bx[i][0])
        parts={nm:split(b[nm],k) for nm in b}
        a_,c_=bx[k]; m_=(a_+c_)/2
        bl=list(bx); bl[k]=(a_,m_); br=list(bx); br[k]=(m_,c_)
        stack.append(({nm:p[0] for nm,p in parts.items()},dep+1,tuple(bl)))
        stack.append(({nm:p[1] for nm,p in parts.items()},dep+1,tuple(br)))
    return dict(certified=(not stack and not deep),closed=closed,open=len(stack),deep=len(deep),sec=round(time.time()-start,1)),deep

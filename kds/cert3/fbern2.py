# float Bernstein with rigorous error tracking (2D)
import numpy as np, flint
from math import comb
from fastbern import rev
def exact_bern(p):
    ctx=p.context(); g=ctx.gens(); one=ctx.from_dict({(0,)*len(g):1})
    r=rev(p,ctx); s=r.compose(*[x+one for x in g],ctx=ctx); d=[int(z) for z in p.degrees()]
    B=np.zeros(tuple(n+1 for n in d),dtype=object)
    for ex,c in s.to_dict().items():
        i,j=[n-int(e) for n,e in zip(d,ex)]
        B[i,j]=flint.fmpq(int(c),comb(d[0],i)*comb(d[1],j))
    M=max(abs(x) for x in B.flat if x!=0)
    A=np.array([[float(x/M) if x!=0 else 0.0 for x in row] for row in B])
    return A   # error per coeff <= 1.2e-16 relative to 1
def split(A,ax):
    A=np.moveaxis(A,ax,0); n=A.shape[0]-1; L=[A[0]]; R=[A[-1]]; T=A.copy()
    for k in range(1,n+1):
        T=0.5*(T[:-1]+T[1:]); L.append(T[0]); R.append(T[-1])
    L=np.array(L); R=np.array(R[::-1])
    return np.moveaxis(L,0,ax),np.moveaxis(R,0,ax)
def fsign(A,err):
    tol=err+1e-300  # absolute: arrays are normalized to max|b|=1 at conversion; fp error <= u(1+sum deg) per item 222
    if A.min()>tol: return 1
    if A.max()<-tol: return -1
    return 0

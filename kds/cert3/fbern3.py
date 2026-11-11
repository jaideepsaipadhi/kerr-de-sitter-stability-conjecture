import numpy as np, flint, itertools
from math import comb
from fastbern import rev
def exact_bern(p,keep=False):
    ctx=p.context(); g=ctx.gens(); one=ctx.from_dict({(0,)*len(g):1})
    r=rev(p,ctx); s=r.compose(*[x+one for x in g],ctx=ctx); d=[int(z) for z in p.degrees()]
    B={}
    for ex,c in s.to_dict().items():
        idx=tuple(n-int(e) for n,e in zip(d,ex))
        den=1
        for n,i in zip(d,idx): den*=comb(n,i)
        B[idx]=flint.fmpq(int(c),den)
    M=max(abs(x) for x in B.values())
    A=np.zeros(tuple(n+1 for n in d))
    for idx,x in B.items(): A[idx]=float(x/M)
    return A*float(M) if keep else A
def split(A,ax):
    A=np.moveaxis(A,ax,0); n=A.shape[0]-1; L=[A[0]]; R=[A[-1]]; T=A.copy()
    for k in range(1,n+1):
        T=0.5*(T[:-1]+T[1:]); L.append(T[0]); R.append(T[-1])
    return np.moveaxis(np.array(L),0,ax),np.moveaxis(np.array(R[::-1]),0,ax)
def fsign(A,err):
    tol=err+1e-300  # absolute: arrays are normalized to max|b|=1 at conversion; fp error <= u(1+sum deg) per item 222
    if A.min()>tol: return 1
    if A.max()<-tol: return -1
    return 0

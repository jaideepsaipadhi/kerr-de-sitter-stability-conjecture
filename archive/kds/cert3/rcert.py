# Recursive exact certificate: p > 0 on (0,1)^3 (open cube), handling edges where p vanishes identically by
# carving epsilon-neighbourhoods and blowing up (two monomial charts, strip the exceptional power), recursively.
import flint, numpy as np, time, sys
from fractions import Fraction
from math import comb, lcm
sys.path.insert(0,'/home/claude/kds/cert2'); from bern_core import split, mn, mx
Q=flint.fmpq_mpoly_ctx.get(('a','b','c'),'lex'); A_,B_,C_=Q.gens(); ONE=Q.from_dict({(0,0,0):1}); G=[A_,B_,C_]
EPS=flint.fmpq(1,8)
LOG=[]
def bern_of(p,box):
    # compose to unit cube on box, then Bernstein (exact ints, positive common scale)
    sub=[ONE*lo+ (hi-lo)*g for g,(lo,hi) in zip(G,box)]
    q=p.compose(*sub,ctx=Q); d=q.to_dict()
    n=[max([e[k] for e in d]+[0]) for k in range(3)]
    A=np.zeros([m+1 for m in n],dtype=object); A.fill(Fraction(0))
    for e,c in d.items(): A[e]=Fraction(int(c.p),int(c.q))
    for k,nn in enumerate(n):
        A=np.moveaxis(A,k,0); Bc=np.empty_like(A)
        for i in range(nn+1):
            acc=A[0]*0
            for j in range(i+1): acc=acc+A[j]*Fraction(comb(i,j),comb(nn,j))
            Bc[i]=acc
        A=np.moveaxis(Bc,0,k)
    L=1
    for x in A.flat: L=lcm(L,x.denominator)
    return np.vectorize(lambda x:int(x*L),otypes=[object])(A)
def vanishing_edges(p):
    out=[]
    for i in range(3):
        for j in range(i+1,3):
            k=3-i-j
            for a in (0,1):
                for b in (0,1):
                    sub=[None]*3; sub[i]=ONE*a; sub[j]=ONE*b; sub[k]=G[k]
                    if p.compose(*sub,ctx=Q)==0: out.append((i,a,j,b))
    return out
def inside(box,nb):
    return all(lo>=nlo and hi<=nhi for (lo,hi),(nlo,nhi) in zip(box,nb))
def strip(q,k):
    d=q.to_dict(); m=min(e[k] for e in d)
    return Q.from_dict({tuple(e[t]-(m if t==k else 0) for t in range(3)):c for e,c in d.items()})
DEAD=[0]
def facestrip(p):
    for k in range(3): p=strip(p,k)
    for g in G:
        h=ONE-g
        while True:
            q,r=divmod(p,h)
            if r==0 and q!=0: p=q
            else: break
    return p
def certify(p,level=0,tmax=50,label=''):
    t0=time.time()
    if DEAD[0] and time.time()>DEAD[0]: LOG.append((label,'global timeout')); return False
    p=facestrip(p)
    edges=vanishing_edges(p) if level<12 else []
    nbs=[]
    for (i,a,j,b) in edges:
        nb=[(flint.fmpq(0),flint.fmpq(1))]*3
        nb[i]=(flint.fmpq(0),EPS) if a==0 else (1-EPS,flint.fmpq(1))
        nb[j]=(flint.fmpq(0),EPS) if b==0 else (1-EPS,flint.fmpq(1))
        nbs.append((nb,(i,a,j,b)))
    verts=[]
    if level<12:
        for a in (0,1):
            for b in (0,1):
                for c in (0,1):
                    if p.compose(ONE*a,ONE*b,ONE*c,ctx=Q)==0 and not any((a,b,c)[i]==ea and (a,b,c)[j]==eb for (i,ea,j,eb) in edges): verts.append((a,b,c))
    vnb=[]
    for vt in verts:
        nb=[((flint.fmpq(0),EPS) if f==0 else (1-EPS,flint.fmpq(1))) for f in vt]; vnb.append((nb,vt))
    allnb=nbs+[(nb,None) for nb,_ in vnb]
    # B&B on the cube minus neighbourhoods
    B0=bern_of(p,((flint.fmpq(0),flint.fmpq(1)),)*3)
    stack=[(B0,((flint.fmpq(0),flint.fmpq(1)),)*3)]; closed=0
    while stack:
        if time.time()-t0>tmax or (DEAD[0] and time.time()>DEAD[0]): LOG.append((label,'timeout')); return False
        Bc,box=stack.pop()
        if any(inside(box,nb) for nb,_ in allnb): continue
        if mn(Bc)>0: closed+=1; continue
        if mx(Bc)<0: LOG.append((label,'NEG',box)); return False
        w=[hi-lo for lo,hi in box]; k=max(range(3),key=lambda t:w[t])
        if w[k]<flint.fmpq(1,2**14): LOG.append((label,'deep',[(float(l),float(h)) for l,h in box])); return False
        lo,hi=box[k]; m=(lo+hi)/2
        L,R=split(Bc,k)
        b1=list(box); b1[k]=(lo,m); b2=list(box); b2[k]=(m,hi); stack+=[(L,tuple(b1)),(R,tuple(b2))]
    # charts for each neighbourhood
    for nb,(i,a,j,b) in nbs:
        k=3-i-j
        xi=[None]*3
        def coord(t,fix): return (ONE*0+EPS*G[t]) if fix==0 else (ONE-EPS*G[t])
        for (u,w_) in ((i,j),(j,i)):     # chart: xi_u = xi_w * alpha  (alpha in slot u)
            sub=[None]*3
            fu=a if u==i else b; fw=b if w_==j else a
            sub[w_]=coord(w_,fw)
            sub[u]=(ONE*0+EPS*G[u]*G[w_]) if fu==0 else (ONE-EPS*G[u]*G[w_])
            sub[k]=G[k]
            q=strip(p.compose(*sub,ctx=Q),w_)
            if not certify(q,level+1,tmax=tmax,label=label+f'/e{(i,a,j,b)}c{u}'): return False
    for nb,vt in vnb:
        for k in range(3):
            sub=[]
            for t in range(3):
                base=(G[t] if t==k else G[t]*G[k])
                sub.append((ONE*0+EPS*base) if vt[t]==0 else (ONE-EPS*base))
            q=strip(p.compose(*sub,ctx=Q),k)
            if not certify(q,level+1,tmax=tmax,label=label+f'/v{vt}c{k}'): return False
    LOG.append((label,'ok',closed)); return True

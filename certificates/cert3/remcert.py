# Bernstein B&B on (0,1)^3 minus the domination-certified neighbourhoods of corners/edges.
import sys,time,numpy as np
from fractions import Fraction as Fr
from math import comb,lcm
sys.path.insert(0,'/home/claude/kds/cert2'); from bern_core import split,mn
def to_bern3(p):
    d=p.degrees(); A=np.empty([k+1 for k in d],dtype=object); A.fill(Fr(0))
    for e,c in p.to_dict().items(): A[tuple(int(i) for i in e)]=Fr(int(c))
    for k,n in enumerate(d):
        A=np.moveaxis(A,k,0); B=np.empty_like(A)
        for i in range(n+1):
            acc=A[0]*0
            for j in range(i+1): acc=acc+A[j]*Fr(comb(i,j),comb(n,j))
            B[i]=acc
        A=np.moveaxis(B,0,k)
    L=1
    for v in A.flat: L=lcm(L,v.denominator)
    return np.vectorize(lambda v:int(v*L),otypes=[object])(A)
def inside(bx,reg):  # reg: list of (lo,hi) per axis
    return all(r[0]<=b[0] and b[1]<=r[1] for b,r in zip(bx,reg))
def run(p,regions,tmax=55):
    t0=time.time(); st=[(to_bern3(p),((Fr(0),Fr(1)),)*3)]; closed=0; skipped=0
    while st:
        if time.time()-t0>tmax: return dict(ok=False,timeout=True,open=len(st),closed=closed)
        A,bx=st.pop()
        if any(inside(bx,r) for r in regions): skipped+=1; continue
        if mn(A)>0: closed+=1; continue
        if min(b[1]-b[0] for b in bx)<Fr(1,2**20): return dict(ok=False,stuck=bx,closed=closed)
        k=max(range(3),key=lambda i:bx[i][1]-bx[i][0]); a,c=bx[k]; m=(a+c)/2
        L,R=split(A,k); bl=list(bx); bl[k]=(a,m); br=list(bx); br[k]=(m,c)
        st.append((L,tuple(bl))); st.append((R,tuple(br)))
    return dict(ok=True,closed=closed,skipped=skipped,sec=round(time.time()-t0,1))
F=Fr
def regions_std(ec=F(1,16),ec2=F(1,64),ee=F(1,256)):
    o,z=F(1),F(0)
    return [((1-ec,o),(z,ec),(1-ec,o)),((1-ec,o),(1-ec,o),(1-ec,o)),((1-ec2,o),(1-ec2,o),(z,ec2)),
            ((z,1-ec),(z,ee),(1-ee,o)),((1-ee,o),(ec,1-ec),(1-ee,o)),((1-ee,o),(1-ee,o),(ec2,1-ec)),
            ((1-ec2,o),(z,ec2),(z,ec2)),((1-ee,o),(z,ee),(ec2,1-ec))]
if __name__=='__main__':
    from ch2 import core3,res
    w,nm,i=sys.argv[1],sys.argv[2],int(sys.argv[3])
    F1,Y,P=res[w]; p=core3(F1 if nm=='F1' else P)[i]
    print(run(p,regions_std()))
def run_resume(p,regions,state,tmax=50):
    import pickle,os
    if os.path.exists(state): st,closed,skipped=pickle.load(open(state,'rb'))
    else: st=[(to_bern3(p),((Fr(0),Fr(1)),)*3)]; closed=skipped=0
    t0=time.time(); minw=Fr(1)
    while st and time.time()-t0<tmax:
        A,bx=st.pop()
        if any(inside(bx,r) for r in regions): skipped+=1; continue
        if mn(A)>0: closed+=1; continue
        w=min(b[1]-b[0] for b in bx); minw=min(minw,w)
        if w<Fr(1,2**20): pickle.dump((st+[(A,bx)],closed,skipped),open(state,'wb')); return dict(ok=False,stuck=bx)
        k=max(range(3),key=lambda i:bx[i][1]-bx[i][0]); a,c=bx[k]; m=(a+c)/2
        L,R=split(A,k); bl=list(bx); bl[k]=(a,m); br=list(bx); br[k]=(m,c)
        st.append((L,tuple(bl))); st.append((R,tuple(br)))
    pickle.dump((st,closed,skipped),open(state,'wb'))
    return dict(ok=not st,open=len(st),closed=closed,skipped=skipped,minw=str(minw))

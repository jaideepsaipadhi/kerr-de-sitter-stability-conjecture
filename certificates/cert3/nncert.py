# Certificate p>0 on the OPEN cube (0,1)^3 by Bernstein subdivision allowing zero coefficients:
# on a box, p = sum b_I B_I with B_I>0 on the relative interior of every face containing index I.
# Accept a box if all b_I>=0 and every face of the box not lying in the cube boundary carries some b_I>0.
import sys,time,itertools,pickle,os,numpy as np
from fractions import Fraction as Fr
from remcert import to_bern3
sys.path.insert(0,'/home/claude/kds/cert2'); from bern_core import split,mn
def ok_box(A,bx):
    if mn(A)<0: return False
    n=A.shape
    for sel in itertools.product((0,1,2),repeat=3):   # 0=lo face,1=hi face,2=free
        if any((s==0 and bx[k][0]==0) or (s==1 and bx[k][1]==1) for k,s in enumerate(sel)): continue
        sl=tuple(0 if s==0 else (n[k]-1 if s==1 else slice(None)) for k,s in enumerate(sel))
        sub=A[sl]
        if not (max(sub.flat) if isinstance(sub,np.ndarray) else sub)>0: return False
    return True
def run(p,state=None,tmax=50,floor=Fr(1,2**16)):
    if state and os.path.exists(state): st,closed=pickle.load(open(state,'rb'))
    else: st=[(to_bern3(p),((Fr(0),Fr(1)),)*3)]; closed=0
    t0=time.time()
    while st and time.time()-t0<tmax:
        A,bx=st.pop()
        if ok_box(A,bx): closed+=1; continue
        if mn(A)<0 and False: pass
        w=min(b[1]-b[0] for b in bx)
        if w<floor: return dict(ok=False,stuck=[str(b) for b in bx],closed=closed)
        k=max(range(3),key=lambda i:bx[i][1]-bx[i][0]); a,c=bx[k]; m=(a+c)/2
        L,R=split(A,k); bl=list(bx); bl[k]=(a,m); br=list(bx); br[k]=(m,c)
        st.append((L,tuple(bl))); st.append((R,tuple(br)))
    if state: pickle.dump((st,closed),open(state,'wb'))
    vol=sum((b[0][1]-b[0][0])*(b[1][1]-b[1][0])*(b[2][1]-b[2][0]) for _,b in st)
    return dict(ok=not st,open=len(st),closed=closed,openvol=float(vol),sec=round(time.time()-t0,1))
if __name__=='__main__':
    from ch2 import core3,res
    w,nm,i=sys.argv[1],sys.argv[2],int(sys.argv[3])
    F1,Y,P=res[w]; p=core3(F1 if nm=='F1' else P)[i]
    print(w,nm,i,run(p,f'nn_{w}_{nm}_{i}.pkl',tmax=int(sys.argv[4])))

import flint, pickle, time, numpy as np, sys
from fractions import Fraction
from math import comb, lcm
sys.path.insert(0,'/home/claude/kds/cert2'); from bern_core import split, mn, mx
t0=time.time()
ctx=flint.fmpz_mpoly_ctx.get(('r0','r2','b','tau'),'lex')
res=pickle.load(open('win_res.pkl','rb'))
def to_arr(f):
    d=f.to_dict(); n2=max(e[1] for e in d); n0=max(e[0] for e in d); n3=max(e[3] for e in d)
    A=np.zeros((n0+1,n2+1,n3+1),dtype=object)
    for e,c in d.items(): A[e[0],n2-e[1],e[3]]+=int(c)   # v^n2 * f(r0,1/v,tau)
    return A
def bern(A):
    A=A.astype(object)
    for k in range(3):
        A=np.moveaxis(A,k,0); n=A.shape[0]-1; Bc=np.empty_like(A)
        for i in range(n+1):
            acc=A[0]*0
            for j in range(i+1): acc=acc+A[j]*Fraction(comb(i,j),comb(n,j))
            Bc[i]=acc
        A=np.moveaxis(Bc,0,k)
    L=1
    for x in A.flat: L=lcm(L,x.denominator)
    return np.vectorize(lambda x:int(x*L),otypes=[object])(A)
def bb(B,tmax=40,maxdep=40):
    st=[(B,0)]; closed=0; start=time.time()
    while st:
        if time.time()-start>tmax: return False,closed,len(st)
        b,d=st.pop()
        if mn(b)>0: closed+=1; continue
        if mx(b)<0 or d>=maxdep: return False,closed,('neg' if mx(b)<0 else 'deep')
        k=int(np.argmax(b.shape)); L,R=split(b,k); st+= [(L,d+1),(R,d+1)]
    return True,closed,0
which=sys.argv[1] if len(sys.argv)>1 else 'Q2'
F1,Y,P=res[which]
for nm,d in (('F1',F1),('P',P)):
    c,fl=ctx.from_dict(d).factor()
    for f,m in fl:
        if sum(f.degrees())<=8:
            continue
        A=to_arr(f); B=bern(A)
        sg=1 if mn(B)>0 else (-1 if mx(B)<0 else 0)
        if sg==0:
            ok,cl,info=bb(B); print(which,nm,'core',f.degrees(),'B&B',ok,cl,info,round(time.time()-t0,1),flush=True)
        else: print(which,nm,'core',f.degrees(),'all Bernstein coeffs sign',sg,round(time.time()-t0,1),flush=True)

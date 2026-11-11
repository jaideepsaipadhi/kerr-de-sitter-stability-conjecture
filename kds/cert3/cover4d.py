import flint,pickle,time,sys,os,json
from fastbern import bern_signs
from fbern3 import exact_bern, split, fsign
RES,NAME,TMAX,DMAX=sys.argv[1],sys.argv[2],float(sys.argv[3]),int(sys.argv[4])
ctx5=flint.fmpz_mpoly_ctx.get(('r0','r2','b','tau','sg'),'lex'); c4=flint.fmpz_mpoly_ctx.get(('x','v','t','s'),'lex')
RES2=pickle.load(open(RES,'rb'))
def tocube(f):
    d=f.to_dict(); n2=max(int(e[1]) for e in d)
    return c4.from_dict({(int(e[0]),n2-int(e[1]),int(e[3]),int(e[4])):int(c) for e,c in d.items()})
def sgnw(p):
    mn,mx,_=bern_signs(p); return 1 if (mn>=0 and mx>0) else (-1 if (mx<=0 and mn<0) else 0)
def prep(d):
    c,fl=ctx5.from_dict(d).factor(); sg=1 if c>0 else -1; cores=[]
    for f,m in fl:
        if m%2==0: continue
        q=tocube(f); w=sgnw(q)
        if w!=0: sg*=w; continue
        cores.append(q)
    return sg,cores
def _scale(q):
    from fastbern import rev
    from math import comb
    ctx=q.context(); g=ctx.gens(); one=ctx.from_dict({(0,)*len(g):1})
    r=rev(q,ctx); s=r.compose(*[x+one for x in g],ctx=ctx); d=[int(z) for z in q.degrees()]
    best=0
    for ex,cc in s.to_dict().items():
        idx=[n-int(e) for n,e in zip(d,ex)]; den=1
        for n,i in zip(d,idx): den*=comb(n,i)
        val=abs(flint.fmpq(int(cc),den))
        if val>best: best=val
    return float(best)
CACHE=f'margin/cache4_{NAME}.pkl'
if os.path.exists(CACHE): S,ITEMS=pickle.load(open(CACHE,'rb'))
else:
    S={}; ITEMS={}
    for k,(s0,facs) in RES2.items():
        its=[]
        for j,(kind,dat) in enumerate(facs):
            if kind=='plain':
                sg,cs=prep(dat); S[(k,j,'v')]=(sg,[exact_bern(q) for q in cs]); its.append(('plain',j))
            else:
                for nm,i in (('F1',0),('Y',1),('P',2)):
                    sg,cs=prep(dat[i]); S[(k,j,nm)]=(sg,[exact_bern(q) for q in cs])
                raw=dat[6]; n2=max(int(e[1]) for e in raw); kmax=max(int(e[2]) for e in raw)
                for kk in range(kmax+1):
                    dk={(int(e[0]),n2-int(e[1]),int(e[3]),int(e[4])):int(cc) for e,cc in raw.items() if int(e[2])==kk}
                    if dk:
                        q=c4.from_dict(dk)
                        from fbern3 import exact_bern as _eb
                        import flint as _fl
                        # recompute true scale: exact_bern normalizes by max |bern coeff|; recover via one coefficient
                        S[(k,j,'c%d'%kk)]=(1,[exact_bern(q,True)])
                its.append(('b',j))
        ITEMS[k]=(s0,its)
    pickle.dump((S,ITEMS),open(CACHE,'wb'))
t0=time.time()
BOX=[None]
def bsg(ps,k,dep):
    s=S[k][0]
    for A in ps[k]:
        g=fsign(A,1e-13*(dep+4))
        if g==0: return 0
        s*=g
    return s
def rule(ps,k,j,dep,F,Yn,Pn):
    f=bsg(ps,(k,j,F),dep); y=bsg(ps,(k,j,Yn),dep)
    if f>0 and y<0: return 1
    if f<0 and y>0: return -1
    p=bsg(ps,(k,j,Pn),dep)
    if p>0 and f!=0: return f
    if p<0 and y!=0: return -y
    return 0
import numpy as np, math
def bplus(x,v):
    r0=x; r2=1/v; r3=-(r0+1+r2); rs=[r0,1,r2,r3]
    E2=sum(rs[i]*rs[j] for i in range(4) for j in range(i+1,4)); E4=r0*r2*r3
    return 2*E4/(E2-math.sqrt(E2*E2-4*E4))
def unred(ps,k,j,dep,box):
    (xl,xh),(vl,vh)=box
    if vl<=0: vl=1e-300
    if xh>=1: xh=1-1e-15
    vals=[bplus(x,v) for x in (xl,xh) for v in (vl,vh)]
    blo=min(vals)*(1-1e-6); bhi=max(vals)*(1+1e-6)
    # rigor of [blo,bhi] assumed via monotonic b+ (TODO certify with q-sign)
    lo=hi=0.0; kk=0
    while (k,j,'c%d'%kk) in S or kk<2:
        if (k,j,'c%d'%kk) in S:
            A=ps[(k,j,'c%d'%kk)][0]; M=np.abs(A).max(); tol=1e-13*(dep+4)*M
            cl,ch=A.min()-tol,A.max()+tol
            pl,ph=blo**kk,bhi**kk
            cand=[cl*pl,cl*ph,ch*pl,ch*ph]; lo+=min(cand); hi+=max(cand)
        kk+=1
        if kk>40: break
    if lo>0: return 1
    if hi<0: return -1
    return 0
def bval(ps,k,j,dep):
    r=rule(ps,k,j,dep,'F1','Y','P')
    if r: return r
    return unred(ps,k,j,dep,BOX[0])
def val(ps,k,dep):
    s0,its=ITEMS[k]; s=s0
    for kind,j in its:
        g=bsg(ps,(k,j,'v'),dep) if kind=='plain' else bval(ps,k,j,dep)
        if g==0: return 0
        s*=g
    return s
def closed(ps,dep):
    for m in (1,2):
        if val(ps,f'R{m}',dep)<0: continue
        if val(ps,f'W1_{m}',dep)>0 and val(ps,f'W2_{m}',dep)>0: continue
        return False
    return True
root={k:v[1] for k,v in S.items()}
def descend(lev,idx):
    cur=root
    for ax in range(4):
        for j in range(lev[ax]-1,-1,-1):
            bit=(idx[ax]>>j)&1; cur={k:[split(A,ax)[bit] for A in v] for k,v in cur.items()}
    return cur
PK=f'margin/c4_{NAME}.json'
todo=[tuple(map(tuple,z)) for z in json.load(open(PK))] if os.path.exists(PK) else [((0,)*4,(0,)*4)]
nc=0; rest=[]
while (todo or rest) and time.time()-t0<TMAX:
    if not todo:
        todo=[r for r in rest if sum(r[0])<DMAX]; rest=[r for r in rest if sum(r[0])>=DMAX]
        if not todo: break
    todo.sort(key=lambda z:-sum(z[0])); lev,idx=todo.pop(); base=sum(lev); st=[(descend(lev,idx),lev,idx)]
    while st:
        if time.time()-t0>TMAX: rest+=[(l,i) for _,l,i in st]; break
        ps,lev,idx=st.pop()
        BOX[0]=[(idx[q]/2**lev[q],(idx[q]+1)/2**lev[q]) for q in (0,1)]
        if closed(ps,sum(lev)): nc+=1; continue
        if sum(lev)>=DMAX or sum(lev)-base>=10: rest.append((lev,idx)); continue
        kk=min(range(4),key=lambda i:lev[i])
        ch=[{},{}]
        for k,v in ps.items():
            sp=[split(A,kk) for A in v]; ch[0][k]=[x[0] for x in sp]; ch[1][k]=[x[1] for x in sp]
        for b in (0,1):
            l2=list(lev); l2[kk]+=1; i2=list(idx); i2[kk]=2*i2[kk]+b
            st.append((ch[b],tuple(l2),tuple(i2)))
rest=todo+rest; json.dump(rest,open(PK,'w'))
print(NAME,'closed',nc,'open',len(rest),'vol',sum(2.0**-sum(l) for l,_ in rest),'maxlev',max([sum(l) for l,_ in rest]+[0]),round(time.time()-t0,1))

import flint,pickle,time,sys,os
from fastbern import bern_signs, half
from fbern2 import exact_bern, split as fsplit, fsign, fsign2
ctx4=flint.fmpz_mpoly_ctx.get(('r0','r2','b','tau'),'lex'); c2=flint.fmpz_mpoly_ctx.get(('x','v'),'lex')
T=pickle.load(open('margin/scg_tau.pkl','rb')); CH=sys.argv[1]; TMAX=float(sys.argv[2]); DMAX=int(sys.argv[3])
X,V=c2.gens(); one=c2.from_dict({(0,0):1}); N=8
def tocube(f):
    d=f.to_dict(); n2=max(int(e[1]) for e in d); q=c2.from_dict({(int(e[0]),n2-int(e[1])):int(c) for e,c in d.items()})
    if CH=='sq':   # x=1-X/8, v=1-V/8 on [0,1]^2
        D=sum(int(z) for z in q.degrees()); out=0*one
        for e,c in q.to_dict().items(): out+=int(c)*(N-X)**int(e[0])*(N-V)**int(e[1])*N**(D-int(e[0])-int(e[1]))
        return out
    D=sum(int(z) for z in q.degrees()); out=0*one
    for e,c in q.to_dict().items():
        a,b=int(e[0]),int(e[1])
        if CH=='EA': out+=int(c)*(N-X)**a*(N-X*V)**b*N**(D-a-b)
        else: out+=int(c)*(N-X*V)**a*(N-X)**b*N**(D-a-b)
    dd=out.to_dict(); k=min(int(e[0]) for e in dd)
    return c2.from_dict({(int(e[0])-k,int(e[1])):c for e,c in dd.items()})
def sgn(p):
    mn,mx,_=bern_signs(p); return 1 if mn>0 else (-1 if mx<0 else 0)
def sgnw(p):
    mn,mx,_=bern_signs(p); return 1 if (mn>=0 and mx>0) else (-1 if (mx<=0 and mn<0) else 0)
S={}
for key,trip in T.items():
    if key[1] not in os.environ.get('CONDS','q0,q1,v1').split(','): continue
    for nm,i in (('F1',0),('Y',1),('P',2)):
        c,fl=ctx4.from_dict(trip[i]).factor(); s=1 if c>0 else -1; cores=[]
        for f,m in fl:
            if m%2==0: continue
            if sum(f.degrees())<=2:
                q=tocube(f); sg=sgnw(q)
                if sg==0: cores.append(q)
                else: s*=sg
            else: cores.append(tocube(f))
        S[key+(nm,)]=(s,cores)
print('built',round(time.time(),1)%1000,flush=True)
def bsg(ps,k):
    s,_=S[k]
    for q in ps[k]:
        g=sgn(q)
        if g==0: return 0
        s*=g
    return s
def pos(ps,tk,cn):
    f=bsg(ps,(tk,cn,'F1'))
    if f>0 and bsg(ps,(tk,cn,'Y'))<0: return True
    p=bsg(ps,(tk,cn,'P'))
    return (p>0 and f>0) or (p<0 and bsg(ps,(tk,cn,'Y'))<0)
def closed(ps):
    for tk in ('Q2','H2','H1'):
        if not (pos(ps,tk,'q0') and pos(ps,tk,'q1')): return False
        if not any(pos(ps,tk,c) for c in ('c2neg','v1','v2','disc') if (tk,c,'F1') in S): return False
    return True
t0=time.time()
FS={k:(v[0],[exact_bern(q) for q in v[1]]) for k,v in S.items()}
print('bern built',round(time.time()-t0,1),flush=True)
def bsg(ps,k,dep):
    s=FS[k][0]
    for A in ps[k]:
        g=fsign2(A,1e-13*(dep+4))
        if g==0: return 0
        s*=g
    return s
def pos(ps,tk,cn,dep):
    f=bsg(ps,(tk,cn,'F1'),dep)
    if f>0 and bsg(ps,(tk,cn,'Y'),dep)<0: return True
    p=bsg(ps,(tk,cn,'P'),dep)
    return (p>0 and f>0) or (p<0 and bsg(ps,(tk,cn,'Y'),dep)<0)
def closed(ps,dep):
    for tk in ('Q2','H2','H1'):
        if not (pos(ps,tk,'q0',dep) and pos(ps,tk,'q1',dep)): return False
        if any(pos(ps,tk,c,dep) for c in ('c2neg','v1','v2','disc') if (tk,c,'F1') in FS): continue
        if all(pos(ps,tk,c,dep) for c in ('mm','pm','mp')): continue
        return False
    return True
import json
import numpy as _np
root={k:[_np.stack([A,_np.abs(A)],-1) for A in v[1]] for k,v in FS.items()}
def descend(lev,idx):
    cur=root
    for ax in (0,1):
        for j in range(lev[ax]-1,-1,-1):
            bit=(idx[ax]>>j)&1
            cur={k:[fsplit(A,ax)[bit] for A in v] for k,v in cur.items()}
    return cur
PK=f'margin/c2f_{CH}.json'
todo=[tuple(map(tuple,z)) for z in json.load(open(PK))] if os.path.exists(PK) else [((0,0),(0,0))]
nc=0; rest=[]
todo.sort(key=lambda z:-sum(z[0]))
while todo and time.time()-t0<TMAX:
    lev,idx=todo.pop(); base=sum(lev); st=[(descend(lev,idx),lev,idx)]
    while st:
        if time.time()-t0>TMAX: rest+=[(l,i) for _,l,i in st]; break
        ps,lev,idx=st.pop()
        if closed(ps,sum(lev)): nc+=1; continue
        if sum(lev)>=DMAX: rest.append((lev,idx)); continue
        kk=0 if lev[0]<=lev[1] else 1
        ch=[{},{}]
        for k,v in ps.items():
            sp=[fsplit(A,kk) for A in v]; ch[0][k]=[x[0] for x in sp]; ch[1][k]=[x[1] for x in sp]
        for b in (0,1):
            l2=list(lev); l2[kk]+=1; i2=list(idx); i2[kk]=2*i2[kk]+b
            st.append((ch[b],tuple(l2),tuple(i2)))
rest=todo+rest
json.dump(rest,open(PK,'w'))
print(CH,'closed',nc,'open',len(rest),'vol',sum(2.0**-sum(l) for l,_ in rest),'maxlev',max([sum(l) for l,_ in rest]+[0]),round(time.time()-t0,1))

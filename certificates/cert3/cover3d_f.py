import flint,pickle,time,sys,os,json
from fastbern import bern_signs
from fbern3 import exact_bern, split, fsign, fsign2
RES,TG,NAME,TMAX,DMAX=sys.argv[1],sys.argv[2].split(','),sys.argv[3],float(sys.argv[4]),int(sys.argv[5])
CH=sys.argv[6] if len(sys.argv)>6 else ''
ctx4=flint.fmpz_mpoly_ctx.get(('r0','r2','b','tau'),'lex'); c3=flint.fmpz_mpoly_ctx.get(('x','v','t'),'lex')
X,V,T3=c3.gens(); one3=c3.from_dict({(0,0,0):1}); N=8
res=pickle.load(open(RES,'rb'))
def tocube(f):
    d=f.to_dict(); n2=max(int(e[1]) for e in d); q=c3.from_dict({(int(e[0]),n2-int(e[1]),int(e[3])):int(c) for e,c in d.items()})
    if not CH: return q
    D=sum(int(z) for z in q.degrees()); out=0*one3
    for e,c in q.to_dict().items():
        a,b_,t_=[int(z) for z in e]
        if CH=='EA': out+=int(c)*(N-X)**a*(N-X*V)**b_*T3**t_*N**(D-a-b_)
        elif CH=='EB': out+=int(c)*(N-X*V)**a*(N-X)**b_*T3**t_*N**(D-a-b_)
        elif CH=='KA': out+=int(c)*(N-X)**a*(X*V)**b_*T3**t_*N**(D-a-b_)
        else: out+=int(c)*(N-X*V)**a*(X)**b_*T3**t_*N**(D-a-b_)
    dd=out.to_dict(); k=min(int(e[0]) for e in dd)
    return c3.from_dict({(int(e[0])-k,int(e[1]),int(e[2])):c for e,c in dd.items()})
def sgnw(p):
    mn,mx,_=bern_signs(p); return 1 if (mn>=0 and mx>0) else (-1 if (mx<=0 and mn<0) else 0)
S={}
for k in TG:
    for nm,i in (('F1',0),('Y',1),('P',2)):
        c,fl=ctx4.from_dict(res[k][i]).factor(); s=1 if c>0 else -1; cores=[]
        for f,m in fl:
            if m%2==0: continue
            q=tocube(f)
            sg=sgnw(q)
            if sg!=0: s*=sg; print('definite core',k,nm,f.degrees(),sg,flush=True); continue
            cores.append(q)
        S[(k,nm)]=(s,[exact_bern(q) for q in cores])
t0=time.time()
def bsg(ps,k,dep):
    s=S[k][0]
    for A in ps[k]:
        g=fsign2(A,1e-13*(dep+4))
        if g==0: return 0
        s*=g
    return s
def closed(ps,dep):
    for k in TG:
        f=bsg(ps,(k,'F1'),dep)
        if f>0 and bsg(ps,(k,'Y'),dep)<0: continue
        p=bsg(ps,(k,'P'),dep)
        if (p>0 and f>0) or (p<0 and bsg(ps,(k,'Y'),dep)<0): continue
        return False
    return True
import numpy as _np
root={k:[_np.stack([A,_np.abs(A)],-1) for A in v[1]] for k,v in S.items()}
def descend(lev,idx):
    cur=root
    for ax in range(3):
        for j in range(lev[ax]-1,-1,-1):
            bit=(idx[ax]>>j)&1; cur={k:[split(A,ax)[bit] for A in v] for k,v in cur.items()}
    return cur
PK=f'margin/c3f_{NAME}{CH}.json'
todo=[tuple(map(tuple,z)) for z in json.load(open(PK))] if os.path.exists(PK) else [((0,0,0),(0,0,0))]
todo=[z for z in todo if CH or not (all(z[1][k]*8>=7*2**z[0][k] for k in (0,1)) or (z[1][0]*8>=7*2**z[0][0] and (z[1][1]+1)*8<=2**z[0][1]))]
todo.sort(key=lambda z:-sum(z[0])); nc=0; rest=[]
while (todo or rest) and time.time()-t0<TMAX:
    if not todo:
        todo=[r for r in rest if sum(r[0])<DMAX]; rest=[r for r in rest if sum(r[0])>=DMAX]
        if not todo: break
        todo.sort(key=lambda z:-sum(z[0]))
    lev,idx=todo.pop(); base=sum(lev); st=[(descend(lev,idx),lev,idx)]
    while st:
        if time.time()-t0>TMAX: rest+=[(l,i) for _,l,i in st]; break
        ps,lev,idx=st.pop()
        if (not CH and (all(idx[k]*8>=7*2**lev[k] for k in (0,1)) or (idx[0]*8>=7*2**lev[0] and (idx[1]+1)*8<=2**lev[1]))) or closed(ps,sum(lev)): nc+=1; continue
        if sum(lev)>=DMAX or sum(lev)-base>=9: rest.append((lev,idx)); continue
        kk=min(range(3),key=lambda i:lev[i])
        ch=[{},{}]
        for k,v in ps.items():
            sp=[split(A,kk) for A in v]; ch[0][k]=[x[0] for x in sp]; ch[1][k]=[x[1] for x in sp]
        for b in (0,1):
            l2=list(lev); l2[kk]+=1; i2=list(idx); i2[kk]=2*i2[kk]+b
            st.append((ch[b],tuple(l2),tuple(i2)))
rest=todo+rest; json.dump(rest,open(PK,'w'))
print(NAME,CH,'closed',nc,'open',len(rest),'vol',sum(2.0**-sum(l) for l,_ in rest),'maxlev',max([sum(l) for l,_ in rest]+[0]),round(time.time()-t0,1))

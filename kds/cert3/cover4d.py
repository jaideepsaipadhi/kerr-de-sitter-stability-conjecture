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
                for nm,i in (('F1',0),('Y',1),('P',2),('F1a',3),('Ya',4),('Pa',5)):
                    sg,cs=prep(dat[i]); S[(k,j,nm)]=(sg,[exact_bern(q) for q in cs])
                its.append(('b',j))
        ITEMS[k]=(s0,its)
    pickle.dump((S,ITEMS),open(CACHE,'wb'))
t0=time.time()
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
def bval(ps,k,j,dep):
    r=rule(ps,k,j,dep,'F1','Y','P')
    if r: return r
    return -rule(ps,k,j,dep,'F1a','Ya','Pa')
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

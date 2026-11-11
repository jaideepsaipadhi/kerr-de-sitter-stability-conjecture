import flint,pickle,time,sys,os,json
from fastbern import bern_signs
from fbern3 import exact_bern, split, fsign
RES,TG,NAME,TMAX,DMAX=sys.argv[1],sys.argv[2].split(','),sys.argv[3],float(sys.argv[4]),int(sys.argv[5])
CH=sys.argv[6] if len(sys.argv)>6 else ''
ctx4=flint.fmpz_mpoly_ctx.get(('r0','r2','b','tau'),'lex'); c3=flint.fmpz_mpoly_ctx.get(('x','v','t'),'lex')
X,V,T3=c3.gens(); one3=c3.from_dict({(0,0,0):1}); N=4 if CH in ('EA','EB') else 8
def tocube(f):
    d=f.to_dict(); n2=max(int(e[1]) for e in d); q=c3.from_dict({(int(e[0]),n2-int(e[1]),int(e[3])):int(c) for e,c in d.items()})
    if not CH: return q
    D=sum(int(z) for z in q.degrees()); out=0*one3
    for e,c in q.to_dict().items():
        a,b_,t_=[int(z) for z in e]
        if CH=='EA': out+=int(c)*(N-X)**a*(N-X*V)**b_*T3**t_*N**(D-a-b_)
        elif CH=='EB': out+=int(c)*(N-X*V)**a*(N-X)**b_*T3**t_*N**(D-a-b_)
        elif CH=='XT1': out+=int(c)*(N-X)**a*(N*V)**b_*(N-X*T3)**t_*N**(D-a-b_-t_)
        elif CH=='XT2': out+=int(c)*(N-X*T3)**a*(N*V)**b_*(N-X)**t_*N**(D-a-b_-t_)
        elif CH=='VT1': out+=int(c)*(N*X)**a*(V)**b_*(N-V*T3)**t_*N**(D-a-b_-t_)*N**0
        elif CH=='VT2': out+=int(c)*(N*X)**a*(V*T3)**b_*(N-V)**t_*N**(D-a-b_-t_)
        elif CH=='KA': out+=int(c)*(N-X)**a*(X*V)**b_*T3**t_*N**(D-a-b_)
        else: out+=int(c)*(N-X*V)**a*(X)**b_*T3**t_*N**(D-a-b_)
    dd=out.to_dict(); sl=1 if CH.startswith('VT') else 0; k=min(int(e[sl]) for e in dd)
    sh=[0,0,0]; sh[sl]=k
    for q in range(3):
        if q!=sl: sh[q]=min(int(e[q]) for e in dd)
    return c3.from_dict({tuple(int(e[q])-sh[q] for q in range(3)):c for e,c in dd.items()})
def sgnw(p):
    mn,mx,_=bern_signs(p); return 1 if (mn>=0 and mx>0) else (-1 if (mx<=0 and mn<0) else 0)
RES2=pickle.load(open(RES,'rb'))
def prep(d):
    c,fl=ctx4.from_dict(d).factor(); sg=1 if c>0 else -1; cores=[]
    for f,m in fl:
        if m%2==0: continue
        q=tocube(f); w=sgnw(q)
        if w!=0: sg*=w; continue
        cores.append(q)
    return sg,cores
CACHE=f'margin/cache_{NAME}{CH}_'+'-'.join(TG)+'.pkl'
if os.path.exists(CACHE): S,ITEMS=pickle.load(open(CACHE,'rb'))
else:
    S={}; ITEMS={}
    for k in TG:
        s0,facs=RES2[k]; its=[]
        for j,(kind,dat) in enumerate(facs):
            if kind=='plain':
                sg,cs=prep(dat); S[(k,j,'v')]=(sg,[exact_bern(q) for q in cs]); its.append(('plain',j))
            else:
                for nm,i in (('F1',0),('Y',1),('P',2)):
                    sg,cs=prep(dat[i]); S[(k,j,nm)]=(sg,[exact_bern(q) for q in cs])
                its.append(('b',j))
        ITEMS[k]=(s0,its)
    
    pickle.dump((S,ITEMS),open(CACHE,'wb'))
ALT=os.environ.get('ALT'); ALTS=[]
if ALT:
    _S,_I,_TG,_R2,NAME0=S,ITEMS,TG,RES2,NAME; SS=dict(S)
    for ai,spec in enumerate(ALT.split(';')):
        altres,alttg=spec.split(':'); TG=alttg.split(','); RES2=pickle.load(open(altres,'rb')); NAME=NAME0+'ALT'+str(ai)+os.path.basename(altres)
        CACHE=f'margin/cache_{NAME}{CH}_'+'-'.join(TG)+'.pkl'
        if os.path.exists(CACHE): S,ITEMS=pickle.load(open(CACHE,'rb'))
        else:
            S={}; ITEMS={}
            for k in TG:
                s0,facs=RES2[k]; its=[]
                for j,(kind,dat) in enumerate(facs):
                    if kind=='plain':
                        sg,cs=prep(dat); S[(k,j,'v')]=(sg,[exact_bern(q) for q in cs]); its.append(('plain',j))
                    else:
                        for nm,i in (('F1',0),('Y',1),('P',2)):
                            sg,cs=prep(dat[i]); S[(k,j,nm)]=(sg,[exact_bern(q) for q in cs])
                        its.append(('b',j))
                ITEMS[k]=(s0,its)
            
            pickle.dump((S,ITEMS),open(CACHE,'wb'))
        
        for k,v in S.items(): SS[('A%d'%ai,)+k]=v
        ALTS.append((TG,ITEMS,('A%d'%ai,)))
    S,ITEMS,TG,RES2,NAME=SS,_I,_TG,_R2,NAME0
t0=time.time()
def bsg(ps,k,dep):
    s=S[k][0]
    for A in ps[k]:
        g=fsign(A,1e-13*(dep+4))
        if g==0: return 0
        s*=g
    return s
def bval(ps,kk,j,dep):
    f=bsg(ps,kk+(j,'F1'),dep); y=bsg(ps,kk+(j,'Y'),dep)
    if f>0 and y<0: return 1
    if f<0 and y>0: return -1
    p=bsg(ps,kk+(j,'P'),dep)
    if p>0 and f!=0: return f
    if p<0 and y!=0: return -y
    return 0
def closed(ps,dep):
    if os.environ.get('PERM2'):
        tgL,itL,pL=ALTS[0]
        return all(closed1(ps,dep,[w],ITEMS,()) or closed1(ps,dep,[l],itL,pL) for w,l in zip(TG,tgL))
    if os.environ.get('PERM'):
        (tgK,itK,pK),(tgL,itL,pL)=ALTS[0],ALTS[1]
        if closed1(ps,dep,tgK,itK,pK): return True
        return all(closed1(ps,dep,[w],ITEMS,()) or closed1(ps,dep,[l],itL,pL) for w,l in zip(TG,tgL))
    if closed1(ps,dep,TG,ITEMS,()): return True
    return any(closed1(ps,dep,tg,it,pre) for tg,it,pre in ALTS)
def closed1(ps,dep,TGx,IT,pre):
    for k in TGx:
        s0,its=IT[k]; s=s0
        for kind,j in its:
            g=bsg(ps,pre+(k,j,'v'),dep) if kind=='plain' else bval(ps,pre+(k,),j,dep)
            if g==0: return False
            s*=g
        if s<=0: return False
    return True
root={k:v[1] for k,v in S.items()}
def descend(lev,idx):
    cur=root
    for ax in range(3):
        for j in range(lev[ax]-1,-1,-1):
            bit=(idx[ax]>>j)&1; cur={k:[split(A,ax)[bit] for A in v] for k,v in cur.items()}
    return cur
PK=f'margin/c3f_{NAME}{CH}.json'
todo=[tuple(map(tuple,z)) for z in json.load(open(PK))] if os.path.exists(PK) else [((0,0,0),(0,0,0))]
todo=todo if os.environ.get('NOEX') else [z for z in todo if not (CH in ('XT1','XT2') and (z[1][1]*4>=3*2**z[0][1] or (z[1][1]+1)*8<=2**z[0][1]))]
todo=todo if os.environ.get('NOEX') else [z for z in todo if CH or not ((z[1][0]*8>=7*2**z[0][0] and z[1][2]*8>=7*2**z[0][2]) or ((z[1][1]+1)*8<=2**z[0][1] and z[1][2]*8>=7*2**z[0][2]) or all(z[1][k]*8>=7*2**z[0][k] for k in (0,1)) or (z[1][0]*8>=7*2**z[0][0] and (z[1][1]+1)*8<=2**z[0][1]))]
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
        if os.environ.get('NOEX') and closed(ps,sum(lev)): nc+=1; continue
        if os.environ.get('NOEX'): pass
        elif (CH in ('XT1','XT2') and (idx[1]*4>=3*2**lev[1] or (idx[1]+1)*8<=2**lev[1])) or (not CH and ((idx[0]*8>=7*2**lev[0] and idx[2]*8>=7*2**lev[2]) or ((idx[1]+1)*8<=2**lev[1] and idx[2]*8>=7*2**lev[2]) or all(idx[k]*8>=7*2**lev[k] for k in (0,1)) or (idx[0]*8>=7*2**lev[0] and (idx[1]+1)*8<=2**lev[1]))) or closed(ps,sum(lev)): nc+=1; continue
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

import flint,pickle,time,sys,os,json
from fastbern import bern_signs
from fbern3 import exact_bern, split, fsign
RES,TG,NAME,TMAX,DMAX=sys.argv[1],sys.argv[2].split(','),sys.argv[3],float(sys.argv[4]),int(sys.argv[5])
CH=sys.argv[6] if len(sys.argv)>6 else ''
ctx4=flint.fmpz_mpoly_ctx.get(('r0','r2','b','tau'),'lex'); c3=flint.fmpz_mpoly_ctx.get(('x','v','t'),'lex')
X,V,T3=c3.gens(); one3=c3.from_dict({(0,0,0):1}); N=4 if CH in ('EA','EB') else 8
def blow02(q,lead,strip=True):
    # q(l,u1,u2): lead='a': l=s, u2=s*w ; lead='b': u2=s, l=s*w ; vars (s,u1,w) in slots (0,1,2)
    out=0*one3
    for e,cc in q.to_dict().items():
        a,b_,t_=[int(z) for z in e]
        if lead=='a': out+=int(cc)*X**(a+t_)*V**b_*T3**t_
        else: out+=int(cc)*(X*T3)**a*V**b_*X**t_
    if not strip: return out
    dd=out.to_dict(); k=min(int(e[0]) for e in dd)
    return c3.from_dict({(int(e[0])-k,int(e[1]),int(e[2])):v for e,v in dd.items()})
def blow12(q,lead,strip=True):
    # q(l,a1,a2): lead 'a': a1=s, a2=s*w ; 'b': a2=s, a1=s*w ; new vars (l,s,w)
    out=0*one3
    for e,cc in q.to_dict().items():
        a,b_,t_=[int(z) for z in e]
        if lead=='a': out+=int(cc)*X**a*V**(b_+t_)*T3**t_
        else: out+=int(cc)*X**a*(V*T3)**b_*V**t_
    if not strip: return out
    dd=out.to_dict(); k=min(int(e[1]) for e in dd)
    return c3.from_dict({(int(e[0]),int(e[1])-k,int(e[2])):v for e,v in dd.items()})
CH0=CH
def tocube(f,n2=None,D=None,strip=True):
    global CH
    if CH0 in ('U3a','U3b') and CH!='U3':
        CH='U3'; q=tocube(f,n2,D,strip); CH=CH0; return blow12(q,CH0[-1],strip)
    if CH0 in ('T3a','T3b') and CH!='T3':
        CH='T3'; q=tocube(f,n2,D,strip); CH=CH0; return blow12(q,CH0[-1],strip)
    if CH0 in ('T1a','T1b') and CH!='T1':
        CH='T1'; q=tocube(f,n2,D,strip); CH=CH0; return blow02(q,CH0[-1],strip)
    d=f.to_dict()
    if n2 is None: n2=max(int(e[1]) for e in d)
    q=c3.from_dict({(int(e[0]),n2-int(e[1]),int(e[3])):int(c) for e,c in d.items()})
    if not CH: return q
    if D is None: D=sum(int(z) for z in q.degrees())
    out=0*one3
    for e,c in q.to_dict().items():
        a,b_,t_=[int(z) for z in e]
        if CH=='EA': out+=int(c)*(N-X)**a*(N-X*V)**b_*T3**t_*N**(D-a-b_)
        elif CH=='EB': out+=int(c)*(N-X*V)**a*(N-X)**b_*T3**t_*N**(D-a-b_)
        elif CH=='U1': out+=int(c)*(N-X)**a*(N-X*V)**b_*(N-X*T3)**t_*N**(D-a-b_-t_)
        elif CH=='U2': out+=int(c)*(N-X*V)**a*(N-X)**b_*(N-X*T3)**t_*N**(D-a-b_-t_)
        elif CH=='U3': out+=int(c)*(N-X*V)**a*(N-X*T3)**b_*(N-X)**t_*N**(D-a-b_-t_)
        elif CH=='T1': out+=int(c)*(N-X)**a*(X*V)**b_*(N-X*T3)**t_*N**(D-a-b_-t_)
        elif CH=='T2': out+=int(c)*(N-X*V)**a*(X)**b_*(N-X*T3)**t_*N**(D-a-b_-t_)
        elif CH=='T3': out+=int(c)*(N-X*V)**a*(X*T3)**b_*(N-X)**t_*N**(D-a-b_-t_)
        elif CH=='XT1': out+=int(c)*(N-X)**a*(N*V)**b_*(N-X*T3)**t_*N**(D-a-b_-t_)
        elif CH=='XT2': out+=int(c)*(N-X*T3)**a*(N*V)**b_*(N-X)**t_*N**(D-a-b_-t_)
        elif CH=='VT1': out+=int(c)*(N*X)**a*(V)**b_*(N-V*T3)**t_*N**(D-a-b_-t_)*N**0
        elif CH=='VT2': out+=int(c)*(N*X)**a*(V*T3)**b_*(N-V)**t_*N**(D-a-b_-t_)
        elif CH=='KA': out+=int(c)*(N-X)**a*(X*V)**b_*T3**t_*N**(D-a-b_)
        else: out+=int(c)*(N-X*V)**a*(X)**b_*T3**t_*N**(D-a-b_)
    if not strip: return out
    dd=out.to_dict(); sl=1 if CH.startswith('VT') else 0; k=min(int(e[sl]) for e in dd)
    sh=[0,0,0]; sh[sl]=k
    for q in range(3):
        if q!=sl: sh[q]=min(int(e[q]) for e in dd)
    return c3.from_dict({tuple(int(e[q])-sh[q] for q in range(3)):c for e,c in dd.items()})
def add_raw(S,k,j,dat):
    if len(dat)<4 or dat[3] is None: return
    raw=dat[3]; n2=max(int(e[1]) for e in raw); kmax=max(int(e[2]) for e in raw)
    fr=ctx4.from_dict(raw)
    D=None
    if CH:
        q0=c3.from_dict({(int(e[0]),n2-int(e[1]),int(e[3])):1 for e in raw}); D=sum(int(z) for z in q0.degrees())
    CK={}
    for kk in range(kmax+1):
        dk={e:cc for e,cc in raw.items() if int(e[2])==kk}
        if dk: CK[kk]=tocube(ctx4.from_dict(dk),n2,D,False)
    if CH:
        sl=1 if CH.startswith('VT') else 0
        if CH0 in ('T3a','T3b','U3a','U3b'):
            for sl_ in (0,1):
                km_=min(min(int(e[sl_]) for e in q.to_dict()) for q in CK.values())
                CK={kk:c3.from_dict({tuple(int(e[w])-(km_ if w==sl_ else 0) for w in range(3)):v for e,v in q.to_dict().items()}) for kk,q in CK.items()}
            sl=2
        km=min(min(int(e[sl]) for e in q.to_dict()) for q in CK.values())
        CK={kk:c3.from_dict({tuple(int(e[w])-(km if w==sl else 0) for w in range(3)):v for e,v in q.to_dict().items()}) for kk,q in CK.items()}
    for kk,q in CK.items(): S[(k,j,'c%d'%kk)]=(1,[exact_bern(q,True)])
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
                add_raw(S,k,j,dat); its.append(('b',j))
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
                        add_raw(S,k,j,dat); its.append(('b',j))
                ITEMS[k]=(s0,its)
            
            pickle.dump((S,ITEMS),open(CACHE,'wb'))
        
        for k,v in S.items(): SS[('A%d'%ai,)+k]=v
        ALTS.append((TG,ITEMS,('A%d'%ai,)))
    S,ITEMS,TG,RES2,NAME=SS,_I,_TG,_R2,NAME0
t0=time.time(); BOX=[None]
BOX=[None]
def bsg(ps,k,dep):
    s=S[k][0]
    for A in ps[k]:
        g=fsign2(A,1e-13*(dep+4))
        if g==0: return 0
        s*=g
    return s
import math
def xv_of(cb):
    # cb: chart coords (c0,c1,c2) -> (x,v)
    a0,a1,a2=cb
    if not CH0: return a0,a1
    if CH0=='U1': return 1-a0/N, 1-a0*a1/N
    if CH0=='U2': return 1-a0*a1/N, 1-a0/N
    if CH0=='U3': return 1-a0*a1/N, 1-a0*a2/N
    if CH0=='U3a': return 1-a0*a1/N, 1-a0*a1*a2/N
    if CH0=='U3b': return 1-a0*a1*a2/N, 1-a0*a1/N
    if CH0=='T1': return 1-a0/N, a0*a1/N
    if CH0 in ('T1a','T1b'):
        lam=a0 if CH0=='T1a' else a0*a2
        return 1-lam/N, lam*a1/N
    if CH0=='T2': return 1-a0*a1/N, a0/N
    if CH0=='T3': return 1-a0*a1/N, a0*a2/N
    if CH0=='T3a': return 1-a0*a1/N, a0*a1*a2/N
    if CH0=='T3b': return 1-a0*a1*a2/N, a0*a1/N
    if CH0=='VT1' or CH0=='VT2': return a0, (a1/N if CH0=='VT1' else a1*a2/N)
    if CH0=='XT1': return 1-a0/N, a1
    if CH0=='XT2': return 1-a0*a2/N, a1
    if CH0=='KA': return 1-a0/N, a0*a1/N
    if CH0=='KB': return 1-a0*a1/N, a0/N
    if CH0=='EA': return 1-a0/N, 1-a0*a1/N
    if CH0=='EB': return 1-a0*a1/N, 1-a0/N
def bplus(x,v):
    r0=x; r2=1/v; r3=-(r0+1+r2); rs=[r0,1,r2,r3]
    E2=sum(rs[i]*rs[jj] for i in range(4) for jj in range(i+1,4)); E4=r0*r2*r3
    return 2*E4/(E2-math.sqrt(E2*E2-4*E4))
import numpy as _np
from mpmath import iv
iv.dps=30
def bplus_iv(x,v):
    # x=r0, v=1/r2 intervals; smooth at v=0
    e4v=-x*(1+v*(1+x)); e2v=v*v*x+v+x*v-(v*(1+x)+1)**2
    Dv=e2v*e2v-4*v*v*e4v
    return 2*e4v/(e2v-iv.sqrt(Dv))
def fsign2(A,err):
    v=A[...,0]; t=err*A[...,1]+1e-300
    if (v>t).all(): return 1
    if (v<-t).all(): return -1
    return 0
M0={}
def unred(ps,kk,j,dep):
    if kk+(j,'c0') not in S and kk+(j,'c1') not in S: return 0
    bx=BOX[0]
    cb=[iv.mpf([lo,hi]) for lo,hi in bx]
    try:
        x,v=xv_of(tuple(cb))
        bb=bplus_iv(x,v); blo=float(bb.a)*(1-1e-12); bhi=float(bb.b)*(1+1e-12)
    except Exception: return 0
    if not (blo>0): return 0
    lo=hi=ab=0.0
    for q in range(60):
        key=kk+(j,'c%d'%q)
        if key not in S: continue
        A=ps[key][0]; tl=1e-13*(dep+4)*A[...,1]
        cl,ch=(A[...,0]-tl).min(),(A[...,0]+tl).max(); pl,ph=blo**q,bhi**q
        cand=[cl*pl,cl*ph,ch*pl,ch*ph]; lo+=min(cand); hi+=max(cand); ab+=max(abs(x) for x in cand)
    if lo>1e-12*ab: return 1
    if hi<-1e-12*ab: return -1
    return 0
def bval(ps,kk,j,dep):
    r=bval0(ps,kk,j,dep)
    return r if r else unred(ps,kk,j,dep)
def bval0(ps,kk,j,dep):
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
root={k:[_np.stack([A,_np.abs(A)],-1) for A in v[1]] for k,v in S.items()}
def descend(lev,idx):
    cur=root
    for ax in range(3):
        for j in range(lev[ax]-1,-1,-1):
            bit=(idx[ax]>>j)&1; cur={k:[split(A,ax)[bit] for A in v] for k,v in cur.items()}
    return cur
PK=f'margin/c3f_{NAME}{CH}.json'
todo=[tuple(map(tuple,z)) for z in json.load(open(PK))] if os.path.exists(PK) else [((0,0,0),(0,0,0))]
todo=todo if os.environ.get('NOEX') else [z for z in todo if not (CH in ('XT1','XT2') and not os.environ.get('XTFULL') and (z[1][1]*8>=7*2**z[0][1] or (z[1][1]+1)*8<=2**z[0][1]))]
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
        BOX[0]=[(idx[q]/2**lev[q],(idx[q]+1)/2**lev[q]) for q in range(3)]
        if os.environ.get('NOEX'):
            ge=lambda k,a,b: idx[k]*b>=a*2**lev[k]; le=lambda k,a,b: (idx[k]+1)*b<=a*2**lev[k]
            if (ge(0,7,8) and le(1,1,8)) or (ge(0,7,8) and ge(2,7,8)) or (le(1,1,8) and ge(2,7,8)) or (ge(0,3,4) and ge(1,3,4)) or closed(ps,sum(lev)): nc+=1; continue
        if os.environ.get('NOEX'): pass
        elif (CH in ('XT1','XT2') and not os.environ.get('XTFULL') and (idx[1]*8>=7*2**lev[1] or (idx[1]+1)*8<=2**lev[1])) or (not CH and ((idx[0]*8>=7*2**lev[0] and idx[2]*8>=7*2**lev[2]) or ((idx[1]+1)*8<=2**lev[1] and idx[2]*8>=7*2**lev[2]) or all(idx[k]*8>=7*2**lev[k] for k in (0,1)) or (idx[0]*8>=7*2**lev[0] and (idx[1]+1)*8<=2**lev[1]))) or closed(ps,sum(lev)): nc+=1; continue
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

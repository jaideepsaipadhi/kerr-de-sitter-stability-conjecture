# generic resumable covering: certify value(target)>0 at b+ for all targets on the open window cube
import flint,pickle,time,sys,os
from fastbern import bern_signs, half
ctx4=flint.fmpz_mpoly_ctx.get(('r0','r2','b','tau'),'lex'); c3=flint.fmpz_mpoly_ctx.get(('x','v','t'),'lex')
RES,TG,NAME,TMAX,DMAX=sys.argv[1],sys.argv[2].split(','),sys.argv[3],float(sys.argv[4]),int(sys.argv[5])
CH=sys.argv[6] if len(sys.argv)>6 else ''
res=pickle.load(open(RES,'rb'))
def tocube(f):
    d=f.to_dict(); n2=max(e[1] for e in d); return c3.from_dict({(e[0],n2-e[1],e[3]):int(c) for e,c in d.items()})
def sgn(p):
    mn,mx,_=bern_signs(p); return 1 if mn>0 else (-1 if mx<0 else 0)
def sgnw(p):
    mn,mx,_=bern_signs(p); return 1 if (mn>=0 and mx>0) else (-1 if (mx<=0 and mn<0) else 0)
SC=f'margin/{NAME}_S.pkl'
S={}
if os.path.exists(SC): S={k:(s0,[c3.from_dict(d) for d in cs]) for k,(s0,cs) in pickle.load(open(SC,'rb')).items()}
for k in (TG if not S else []):
    for nm,i in (('F1',0),('Y',1),('P',2)):
        c,fl=ctx4.from_dict(res[k][i]).factor(); s=1 if c>0 else -1; cores=[]
        for f,m in fl:
            if m%2==0: continue
            if sum(f.degrees())<=3:
                sg=sgnw(tocube(f)); assert sg!=0,(k,nm,f); s*=sg
            else: cores.append(tocube(f))
        S[(k,nm)]=(s,cores)
pickle.dump({k:(s0,[q.to_dict() for q in cs]) for k,(s0,cs) in S.items()},open(SC,'wb'))
def boxsign(cs,s):
    for q in cs:
        g=sgn(q)
        if g==0: return 0
        s*=g
    return s
def closed(ps):
    bs=lambda k,nm: boxsign(ps[(k,nm)],S[(k,nm)][0])
    for k in TG:
        f=bs(k,'F1')
        if f>0 and bs(k,'Y')<0: continue
        p=bs(k,'P')
        if (p>0 and f>0) or (p<0 and bs(k,'Y')<0): continue
        return False
    return True
import json
PK=f'margin/{NAME}32{CH.replace(":","_").replace(",","-")}.json'; t0=time.time()
CH=sys.argv[6] if len(sys.argv)>6 else ''
X,V,T3=c3.gens(); one3=c3.from_dict({(0,0,0):1})
def chart2(q,sub):
    # q(sigma,beta,tau) on chart; sub '1': sigma=s/8, tau=s*g/8 ; '2': tau=s/8, sigma=s*g/8 ; vars (s,beta,g)=(X,V,T3)
    N=8; D=sum(int(z) for z in q.degrees()); out=0*one3
    for e,c in q.to_dict().items():
        e=[int(z) for z in e]
        if sub=='1': out+=int(c)*X**e[0]*V**e[1]*(X*T3)**e[2]*N**(D-e[0]-e[2])
        else: out+=int(c)*(X*T3)**e[0]*V**e[1]*X**e[2]*N**(D-e[0]-e[2])
    dd=out.to_dict(); k=min(int(e[0]) for e in dd)
    return c3.from_dict({(int(e[0])-k,int(e[1]),int(e[2])):c for e,c in dd.items()})
def blow(q,i,j,lead):
    # blow up (var_i,var_j)->0 with scale 1/8: lead var = s/8, other = s*g/8 ; s stored in slot i if lead==i else slot j
    N=8; G=list(c3.gens()); D=sum(int(z) for z in q.degrees()); out=0*one3
    sl=lead; ot=j if lead==i else i
    for e,c in q.to_dict().items():
        e=[int(z) for z in e]; term=int(c)
        for k in range(3):
            if k==sl: term=term*G[sl]**e[k]
            elif k==ot: term=term*(G[sl]*G[ot])**e[k]
            else: term=term*G[k]**e[k]
        out+=term*N**(D-e[sl]-e[ot])
    dd=out.to_dict(); k=min(int(e[sl]) for e in dd)
    return c3.from_dict({tuple(int(e[t])-(k if t==sl else 0) for t in range(3)):c for e,c in dd.items()})
def chart(q):
    if not CH: return q
    if ':' in CH:
        base,spec=CH.split(':'); q=chart(q) if False else chart_base(q,base[:2])
        if len(base)==3: q=chart2(q,base[2])
        for tok in spec.split(','):
            i,j,l=int(tok[0]),int(tok[1]),int(tok[2]); q=blow(q,i,j,l)
        return q
    if len(CH)==3:
        return chart2(chart_base(q,CH[:2]),CH[2])
    return chart_base(q,CH)
def chart_base(q,CH):
    N=32
    if CH=='EA': Xs,Vs=(N-X),(N-X*V)     # x=1-sigma/N, v=1-sigma*beta/N ; (sigma,beta,tau)=(X,V,T3)
    else: Xs,Vs=(N-X*V),(N-X)
    D=sum(int(z) for z in q.degrees()); out=0*one3
    for e,c in q.to_dict().items():
        e=[int(z) for z in e]; out+=int(c)*Xs**e[0]*Vs**e[1]*T3**e[2]*N**(e[2]+D-e[0]-e[1]-e[2])
    dd=out.to_dict(); k=min(int(e[0]) for e in dd)
    return c3.from_dict({(int(e[0])-k,int(e[1]),int(e[2])):c for e,c in dd.items()})
root={k:[chart(q) for q in v[1]] for k,v in S.items()}
def excluded(lev,idx):
    if CH in ('EA','EB'): return (idx[0]*8+8<=2**lev[0] and idx[2]*8+8<=2**lev[2]) or (idx[0]*8+8<=2**lev[0] and idx[1]*8+8<=2**lev[1])
    if CH: return False
    return all(idx[k]*32>=31*2**lev[k] for k in (0,1))
def descend(lev,idx):
    out={k:list(v) for k,v in root.items()}
    for ax in range(3):
        for j in range(lev[ax]-1,-1,-1):
            bit=(idx[ax]>>j)&1; out={k:[half(q,ax,bit) for q in v] for k,v in out.items()}
    return out
todo=[tuple(map(tuple,x)) for x in json.load(open(PK))] if os.path.exists(PK) else [((0,0,0),(0,0,0))]
todo=[z for z in todo if not excluded(*z)]
todo.sort(key=lambda z:-sum(z[0])); nc=0; rest=[]
while (todo or rest) and time.time()-t0<TMAX:
    if not todo:
        todo=[r for r in rest if sum(r[0])<DMAX]; keep=[r for r in rest if sum(r[0])>=DMAX]; rest=keep
        if not todo: break
        todo.sort(key=lambda z:-sum(z[0]))
    lev,idx=todo.pop(); base=sum(lev); st=[(descend(lev,idx),lev,idx)]
    while st:
        if time.time()-t0>TMAX: rest+=[(l,i) for _,l,i in st]; break
        st.sort(key=lambda z:-sum(z[1])); ps,lev,idx=st.pop()
        if sum(lev)-base>6: rest.append((lev,idx)); continue
        if excluded(lev,idx) or closed(ps): nc+=1; continue
        if sum(lev)>=DMAX: rest.append((lev,idx)); continue
        kk=min(range(3),key=lambda i:lev[i])
        for b in (0,1):
            l2=list(lev); l2[kk]+=1; i2=list(idx); i2[kk]=2*i2[kk]+b
            st.append(({k:[half(q,kk,b) for q in v] for k,v in ps.items()},tuple(l2),tuple(i2)))
rest=todo+rest
json.dump(rest,open(PK+'.tmp','w')); os.replace(PK+'.tmp',PK)
print(NAME,'closed',nc,'open',len(rest),'vol',sum(2.0**-sum(l) for l,_ in rest),'maxlev',max((sum(l) for l,_ in rest),default=0),round(time.time()-t0,1))

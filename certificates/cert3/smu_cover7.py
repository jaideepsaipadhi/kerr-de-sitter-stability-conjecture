import flint,pickle,time,sys
from fastbern import bern_signs, half
ctx4=flint.fmpz_mpoly_ctx.get(('r0','r2','b','tau'),'lex')
c3=flint.fmpz_mpoly_ctx.get(('x','v','t'),'lex')
res=pickle.load(open('margin/smu_res.pkl','rb'))
def tocube(f):
    d=f.to_dict(); n2=max(e[1] for e in d); return c3.from_dict({(e[0],n2-e[1],e[3]):int(c) for e,c in d.items()})
def sgn(p):
    mn,mx,_=bern_signs(p); return 1 if mn>0 else (-1 if mx<0 else 0)
def sgnw(p):  # weak: >=0 with some >0
    mn,mx,_=bern_signs(p); return 1 if (mn>=0 and mx>0) else (-1 if (mx<=0 and mn<0) else 0)
# for each (target, F1/P): constant sign s0 from definite small factors, and list of core factors (odd mult)
S={}
for k in ('G','Q2','H2','H1'):
    for nm,i in (('F1',0),('Y',1),('P',2)):
        c,fl=ctx4.from_dict(res[k][i]).factor(); s=1 if c>0 else -1; cores=[]
        for f,m in fl:
            if m%2==0: continue
            q=tocube(f)
            if sum(f.degrees())<=3:
                sg=sgnw(q); assert sg!=0,(k,nm,f); s*=sg
            else: cores.append(q)
        S[(k,nm)]=(s,cores)
print({k:(v[0],[p.degrees() for p in v[1]]) for k,v in S.items()},flush=True)
keys=list(S); 
def boxsign(cs,s):
    for q in cs:
        g=sgn(q)
        if g==0: return 0
        s*=g
    return s

import json,os
def descend(polys,lev,idx):
    out={k:list(v) for k,v in polys.items()}
    for ax in range(3):
        for j in range(lev[ax]-1,-1,-1):
            bit=(idx[ax]>>j)&1
            out={k:[half(q,ax,bit) for q in v] for k,v in out.items()}
    return out
root={k:S[k][1] for k in keys}
TMAX=float(sys.argv[1]); DMAX=int(sys.argv[2]); FN='margin/smu_open.json'
PK='margin/smu_state.pkl'
if os.path.exists(PK):
    todo=[(None if pd is None else {k:[c3.from_dict(d) for d in v] for k,v in pd.items()},tuple(l),tuple(i)) for pd,l,i in pickle.load(open(PK,'rb'))]
else:
    todo=[(None,tuple(l),tuple(i)) for l,i in json.load(open(FN))]
t0=time.time(); st=[]; nG=nL=0; openb=[]; leftover=[]
while todo and time.time()-t0<TMAX:
    todo.sort(key=lambda z:-sum(z[1])); it=todo.pop(); st=[(it[0] if it[0] is not None else descend(root,it[1],it[2]),it[1],it[2])]
    while st:
        if time.time()-t0>TMAX: leftover+=st; st=[]; break
        st.sort(key=lambda z:-sum(z[1])); ps,lev2,idx2=st.pop()
        bs=lambda k,nm: boxsign(ps[(k,nm)],S[(k,nm)][0])
        gF=bs('G','F1'); gY=bs('G','Y')
        if gF<0 and gY>0: nG+=1; continue
        gP=bs('G','P')
        if (gP>0 and gF<0) or (gP<0 and gY>0): nG+=1; continue
        def lpos(k):
            f=bs(k,'F1')
            if f>0 and bs(k,'Y')<0: return True
            p=bs(k,'P')
            if p>0: return f>0
            if p<0: return bs(k,'Y')<0
            return False
        if all(lpos(k) for k in ('Q2','H2','H1')): nL+=1; continue
        if sum(lev2)>=DMAX: openb.append((ps,lev2,idx2)); continue
        kk=min(range(3),key=lambda i:lev2[i])
        for r in (0,1):
            l2=list(lev2); l2[kk]+=1; i2=list(idx2); i2[kk]=2*i2[kk]+r
            st.append(({k:[half(q,kk,r) for q in v] for k,v in ps.items()},tuple(l2),tuple(i2)))
rest=todo+leftover+openb
pickle.dump([(None if pd is None else {k:[q.to_dict() for q in v] for k,v in pd.items()},l,i) for pd,l,i in rest],open(PK+'.tmp','wb')); os.replace(PK+'.tmp',PK)
vol=sum(2.0**-sum(l) for _,l,_ in rest)
print('closed G',nG,'L',nL,'remaining',len(rest),'vol',vol,'maxlev',max((sum(l) for _,l,_ in rest),default=0),round(time.time()-t0,1))

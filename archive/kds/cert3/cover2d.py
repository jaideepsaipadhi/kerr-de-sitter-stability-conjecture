import flint,pickle,time,sys,os
from fastbern import bern_signs, half
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
t0=time.time(); st=[({k:v[1] for k,v in S.items()},(0,0),(0,0))]; nc=0; op=[]
while st and time.time()-t0<TMAX:
    st.sort(key=lambda z:-sum(z[1])); ps,lev,idx=st.pop()
    if closed(ps): nc+=1; continue
    if sum(lev)>=DMAX: op.append((lev,idx)); continue
    kk=0 if lev[0]<=lev[1] else 1
    for b in (0,1):
        l2=list(lev); l2[kk]+=1; i2=list(idx); i2[kk]=2*i2[kk]+b
        st.append(({k:[half(q,kk,b) for q in v] for k,v in ps.items()},tuple(l2),tuple(i2)))
rest=[(l,i) for _,l,i in st]+op
print(CH,'closed',nc,'open',len(rest),'vol',sum(2.0**-sum(l) for l,_ in rest),round(time.time()-t0,1))
print(sorted([(round(i[0]/2**l[0],3),round((i[0]+1)/2**l[0],3),round(i[1]/2**l[1],3),round((i[1]+1)/2**l[1],3)) for l,i in rest])[:40])

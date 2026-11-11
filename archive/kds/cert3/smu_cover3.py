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
        if k!='G' and nm=='Y': continue
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
t0=time.time(); st=[({k:S[k][1] for k in keys},(0,0,0),(0,0,0))]; nG=nL=0; openb=[]
TMAX=float(sys.argv[1]); DMAX=int(sys.argv[2])
while st and time.time()-t0<TMAX:
    ps,lev,idx=st.pop()
    gP=boxsign(ps[('G','P')],S[('G','P')][0])
    if gP>0 and boxsign(ps[('G','F1')],S[('G','F1')][0])<0: nG+=1; continue
    if gP<0 and boxsign(ps[('G','Y')],S[('G','Y')][0])>0: nG+=1; continue
    if all(boxsign(ps[(k,nm)],S[(k,nm)][0])>0 for k in ('Q2','H2','H1') for nm in ('F1','P')): nL+=1; continue
    if sum(lev)>=DMAX: openb.append((lev,idx)); continue
    kk=min(range(3),key=lambda i:lev[i])
    for r in (0,1):
        l2=list(lev); l2[kk]+=1; i2=list(idx); i2[kk]=2*i2[kk]+r
        st.append(({k:[half(q,kk,r) for q in v] for k,v in ps.items()},tuple(l2),tuple(i2)))
vol=lambda L: 2.0**-sum(L)
print('closed G',nG,'closed L',nL,'open(depth)',len(openb),'pending',len(st),'openvol',sum(vol(l) for l,_ in openb)+sum(vol(l) for _,l,_ in st),round(time.time()-t0,1))
print('open sample',openb[:10])

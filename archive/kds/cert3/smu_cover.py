import flint,pickle,time,sys
from fastbern import bern_signs, half
ctx4=flint.fmpz_mpoly_ctx.get(('r0','r2','b','tau'),'lex')
c3=flint.fmpz_mpoly_ctx.get(('x','v','t'),'lex')
res=pickle.load(open('margin/smu_res.pkl','rb'))
def tocube(d):
    n2=max(e[1] for e in d); return c3.from_dict({(e[0],n2-e[1],e[3]):int(c) for e,c in d.items()})
polys={k:(tocube(v[0]),tocube(v[2])) for k,v in res.items()}
def sgn(p):
    mn,mx,_=bern_signs(p); return 1 if mn>0 else (-1 if mx<0 else 0)
t0=time.time(); st=[(polys,(0,0,0),(0,0,0))]; nG=nL=0; openb=[]
TMAX=float(sys.argv[1]); DMAX=int(sys.argv[2])
while st and time.time()-t0<TMAX:
    ps,lev,idx=st.pop()
    if sgn(ps['G'][0])<0 and sgn(ps['G'][1])>0: nG+=1; continue
    if all(sgn(ps[k][0])>0 and sgn(ps[k][1])>0 for k in ('Q2','H2','H1')): nL+=1; continue
    if sum(lev)>=DMAX: openb.append((lev,idx)); continue
    k=min(range(3),key=lambda i:lev[i])
    for r in (0,1):
        l2=list(lev); l2[k]+=1; i2=list(idx); i2[k]=2*i2[k]+r
        st.append(({kk:(half(a,k,r),half(b,k,r)) for kk,(a,b) in ps.items()},tuple(l2),tuple(i2)))
vol=lambda L: 2.0**-sum(L)
print('closed G',nG,'closed L',nL,'open(depth)',len(openb),'pending',len(st),'openvol',sum(vol(l) for l,_ in openb)+sum(vol(l) for _,l,_ in st),round(time.time()-t0,1))
pickle.dump((st,openb),open('margin/smu_state.pkl','wb'))

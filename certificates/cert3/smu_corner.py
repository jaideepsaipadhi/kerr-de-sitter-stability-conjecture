import flint,pickle,time,sys
exec(open('smu_cover7.py').read().split("keys=list(S)")[0])
keys=list(S)
cb=flint.fmpz_mpoly_ctx.get(('s','be','ga'),'lex'); s,be,ga=cb.gens(); one=cb.from_dict({(0,0,0):1})
N=256  # s = sigma/N
def chart(q,ch):
    # substitute into c3 poly q(x,v,t), with s=sigma/N : multiply through by N^deg
    if ch=='A': X,V,T=(N-s),s*be,(N-s*ga)
    elif ch=='B': X,V,T=(N-s*be),s,(N-s*ga)
    else: X,V,T=(N-s*be),s*ga,(N-s)
    d=q.to_dict(); out=0*one
    D=sum(int(x) for x in q.degrees())
    for e,c in d.items():
        e=[int(z) for z in e]
        out+=int(c)*X**e[0]*V**e[1]*T**e[2]*N**(D-e[0]-e[1]-e[2])
    dd=out.to_dict(); k=min(int(e[0]) for e in dd)
    return cb.from_dict({(int(e[0])-k,int(e[1]),int(e[2])):c for e,c in dd.items()}),k
import os
ch=sys.argv[1]; t0=time.time()
C=None if os.path.exists(f'margin/corner_{ch}.pkl') else {k:[chart(q,ch)[0] for q in v[1]] for k,v in S.items()}
print('built',ch,round(time.time()-t0,1),flush=True)
def boxsign(cs,s):
    for q in cs:
        g=sgn(q)
        if g==0: return 0
        s*=g
    return s
def bs(ps,k,nm): return boxsign(ps[(k,nm)],S[(k,nm)][0])
def closed(ps):
    gF=bs(ps,'G','F1'); gY=bs(ps,'G','Y')
    if gF<0 and gY>0: return 'G'
    gP=bs(ps,'G','P')
    if (gP>0 and gF<0) or (gP<0 and gY>0): return 'G'
    def lpos(k):
        f=bs(ps,k,'F1')
        if f>0 and bs(ps,k,'Y')<0: return True
        p=bs(ps,k,'P')
        return (p>0 and f>0) or (p<0 and bs(ps,k,'Y')<0)
    return 'L' if all(lpos(k) for k in ('Q2','H2','H1')) else None
import os
PK=f'margin/corner_{ch}.pkl'
if os.path.exists(PK): st=[({k:[cb.from_dict(d) for d in v] for k,v in pd.items()},tuple(l)) for pd,l in pickle.load(open(PK,'rb'))]
else: st=[(C,(0,0,0))]
nc=0; op=[]
while st and time.time()-t0<float(sys.argv[2]):
    st.sort(key=lambda z:-sum(z[1])); ps,lev=st.pop()
    r=closed(ps)
    if r: nc+=1; continue
    if sum(lev)>=int(sys.argv[3]): op.append((ps,lev)); continue
    kk=min(range(3),key=lambda i:lev[i])
    for b in (0,1):
        l2=list(lev); l2[kk]+=1
        st.append(({k:[half(q,kk,b) for q in v] for k,v in ps.items()},tuple(l2)))
rest=st+op
pickle.dump([({k:[q.to_dict() for q in v] for k,v in pd.items()},l) for pd,l in rest],open(PK+'.tmp','wb')); os.replace(PK+'.tmp',PK)
print(ch,'closed',nc,'open',len(rest),'vol',sum(2.0**-sum(l) for _,l in rest),'maxlev',max((sum(l) for _,l in rest),default=0),round(time.time()-t0,1))

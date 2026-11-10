import flint, pickle, numpy as np, time, sys
t0=time.time()
ctx4=flint.fmpz_mpoly_ctx.get(('r0','r2','b','tau'),'lex')
c3=flint.fmpz_mpoly_ctx.get(('x','v','t'),'lex'); X,Vv,T=c3.gens()
res=pickle.load(open('win_res.pkl','rb'))
def core3(d):
    c,fl=ctx4.from_dict(d).factor(); out=[]
    for f,m in fl:
        if sum(f.degrees())<=8: continue
        dd=f.to_dict(); n2=max(e[1] for e in dd)
        out.append(c3.from_dict({(e[0],n2-e[1],e[3]):int(c) for e,c in dd.items()}))
    return out
cA=flint.fmpz_mpoly_ctx.get(('x','mu','nu'),'lex'); x,mu,nu=cA.gens()
def chart(p,which):
    if which=='A': q=p.compose(x, mu*nu, 1-nu, ctx=cA)
    else: q=p.compose(x, nu, 1-mu*nu, ctx=cA)
    d=q.to_dict(); k=min(e[2] for e in d); return cA.from_dict({(e[0],e[1],e[2]-k):c for e,c in d.items()}),k
def fe(p):
    d=p.to_dict(); E=np.array(list(d.keys()),float); C=np.array([float(int(c)) for c in d.values()])
    return lambda a,b_,c: (C*a**E[:,0]*b_**E[:,1]*c**E[:,2]).sum()
if __name__=='__main__':
    which=sys.argv[1]; nm=sys.argv[2]
    F1,Y,P=res[which]; d=F1 if nm=='F1' else P
    for p in core3(d):
        for ch in 'AB':
            q,k=chart(p,ch); f=fe(q); g=np.linspace(0,1,9)
            vals=np.array([[[f(a,b_,c) for c in g] for b_ in g] for a in g])
            z=np.argwhere(np.abs(vals)<=1e-12*np.abs(vals).max())
            print(which,nm,ch,'k',k,'min sign',np.sign(vals.min()),'zeros at (x,mu,nu) grid:',[tuple(g[i] for i in zz) for zz in z][:8],round(time.time()-t0,1),flush=True)

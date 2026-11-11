import flint, pickle, time, numpy as np
t0=time.time()
ctx=flint.fmpz_mpoly_ctx.get(('r0','r2','b','tau'),'lex')
res=pickle.load(open('win_res.pkl','rb'))
def fe(d):
    E=np.array(list(d.keys()),dtype=float); C=np.array([float(int(c)) for c in d.values()])
    def f(x,y,t): return (C[None,:]*x[:,None]**E[:,0]*y[:,None]**E[:,1]*t[:,None]**E[:,3]).sum(1)
    return f
rng=np.random.default_rng(0); X=rng.uniform(0.02,0.98,3000); V=rng.uniform(0.02,0.98,3000); T=rng.uniform(0,1,3000)
for k,(F1,Y,P) in res.items():
    s=[np.sign(fe(d)(X,1/V,T)) for d in (F1,P,Y)]
    print(k,'F1',set(s[0]),'P',set(s[1]),'Y',set(s[2]),flush=True)
    p=ctx.from_dict(P); c,fl=p.factor()
    print('  P factors (deg,mult):',[(f.degrees(),m) for f,m in fl],'const',c, round(time.time()-t0,1),flush=True)
    q=ctx.from_dict(F1); c,fl=q.factor(); print('  F1 factors:',[(f.degrees(),m) for f,m in fl],'const',c,flush=True)

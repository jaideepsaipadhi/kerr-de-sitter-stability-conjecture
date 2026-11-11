import flint,pickle,numpy as np
ctx5=flint.fmpz_mpoly_ctx.get(('r0','r2','b','tau','sg'),'lex')
R=pickle.load(open('margin/sbe2_res.pkl','rb'))
s0,f=R['W1_2']; kind,dat=f[0]
for nm,i in (('F1',0),('Y',1),('P',2)):
    p=ctx5.from_dict(dat[i]); c,fl=p.factor()
    for g,m in fl:
        if m%2==0 or sum(g.degrees())<=2: continue
        d=g.to_dict(); n2=max(int(e[1]) for e in d)
        top={(int(e[0]),int(e[3]),int(e[4])):int(cc) for e,cc in d.items() if int(e[1])==n2}
        X=np.linspace(0.01,0.99,15); vals=np.array([sum(cc*x**e[0]*t**e[1]*s**e[2] for e,cc in top.items()) for x in X for t in (0.2,0.8) for s in (0.1,0.5,0.9)],dtype=float)
        print(nm,g.degrees(),'v=0 face sign',np.sign(vals.min()),np.sign(vals.max()),'zeros',(np.abs(vals)<1e-12*np.abs(vals).max()).sum())

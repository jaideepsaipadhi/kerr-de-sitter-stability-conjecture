import sys; sys.argv=['x','margin/sal2_res.pkl','Q2,H2,H1','sal2','0','40','EA']
exec(open('cover3d_g.py').read().split("t0=time.time()")[0].replace("S[(k,j,'v')]=(sg,[exact_bern(q) for q in cs])","S[(k,j,'v')]=(sg,cs)").replace("S[(k,j,nm)]=(sg,[exact_bern(q) for q in cs])","S[(k,j,nm)]=(sg,cs)"))
import numpy as np
for k in ('H2','H1'):
  s0,its=ITEMS[k]
  if 1:
    Vs=np.linspace(0.0,1.0,101)
    for kind,j in its:
        for nm in (('v',) if kind=='plain' else ('F1','Y','P')):
            sg,cs=S[(k,j,nm)]
            for q in cs:
                f0={(int(e[1]),int(e[2])):int(c) for e,c in q.to_dict().items() if int(e[0])==0}
                vals=np.array([[sum(c*v**e[0]*tt**e[1] for e,c in f0.items()) for v in Vs] for tt in (0.05,0.5,0.95)],dtype=float)
                ch=[[round(Vs[i],3) for i in range(100) if np.sign(r[i])!=np.sign(r[i+1])] for r in vals]
                print(j,nm,q.degrees(),'sign changes in beta (tau=.05,.5,.95)',ch)

import sys; sys.argv=['x','margin/sal2_res.pkl','Q2,H2,H1','sal2','0','40','XT1']
exec(open('cover3d_g.py').read().split("t0=time.time()")[0].replace("S[(k,j,'v')]=(sg,[exact_bern(q) for q in cs])","S[(k,j,'v')]=(sg,cs)").replace("S[(k,j,nm)]=(sg,[exact_bern(q) for q in cs])","S[(k,j,nm)]=(sg,cs)"))
import numpy as np
k='Q2'; s0,its=ITEMS[k]
Vs=np.linspace(0.70,0.875,71)
for kind,j in its:
    for nm in (('v',) if kind=='plain' else ('F1','Y','P')):
        sg,cs=S[(k,j,nm)]
        for q in cs:
            f0={(int(e[1]),int(e[2])):int(c) for e,c in q.to_dict().items() if int(e[0])==0}
            vals=np.array([sum(c*v**e[0]*0.5**e[1] for e,c in f0.items()) for v in Vs],dtype=float)
            ch=[round(Vs[i],4) for i in range(70) if np.sign(vals[i])!=np.sign(vals[i+1])]
            print(j,nm,q.degrees(),'sign changes at v',ch)

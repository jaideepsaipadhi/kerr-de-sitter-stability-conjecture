import sys; sys.argv=['x','margin/sal2_res.pkl','Q2,H2,H1','sal2','0','40','XT1']
exec(open('cover3d_g.py').read().split("t0=time.time()")[0].replace("S[(k,j,'v')]=(sg,[exact_bern(q) for q in cs])","S[(k,j,'v')]=(sg,cs)").replace("S[(k,j,nm)]=(sg,[exact_bern(q) for q in cs])","S[(k,j,nm)]=(sg,cs)"))
import numpy as np
k='Q2'; s0,its=ITEMS[k]
for kind,j in its:
    if kind!='b': continue
    for nm in ('F1','Y','P'):
        sg,cs=S[(k,j,nm)]
        for q in cs:
            f0={(int(e[1]),int(e[2])):int(c) for e,c in q.to_dict().items() if int(e[0])==0}
            vals=[sum(c*v**e[0]*g**e[1] for e,c in f0.items()) for v in np.linspace(0.70,0.875,36) for g in (0.2,0.5,0.8)]
            vv=np.array(vals,dtype=float); print(j,nm,q.degrees(),'face sign range',np.sign(vv.min()),np.sign(vv.max()))

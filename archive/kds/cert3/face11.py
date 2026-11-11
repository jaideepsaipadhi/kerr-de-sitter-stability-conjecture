import sys,os; sys.argv=['x','margin/sal2_res.pkl','H2,H1','sal2','0','40','EB']
src=open('cover3d_h.py').read().split("t0=time.time()")[0]
src=src.replace("[exact_bern(q) for q in cs]","cs").replace("pickle.dump((S,ITEMS),open(CACHE,'wb'))","pass").replace("CACHE=f'margin/cache_","CACHE=f'margin/nocache_")
exec(src)
import numpy as np
for k in ('H2',):
    s0,its=ITEMS[k]
    for kind,j in its:
        if kind!='b': continue
        for nm in ('F1','Y','P'):
            for q in S[(k,j,nm)][1]:
                f0={(int(e[0]),int(e[2])):int(c) for e,c in q.to_dict().items() if int(e[1])==0}
                if not f0: print(j,nm,'face identically 0'); continue
                G=np.array([[sum(c*s**e[0]*t**e[1] for e,c in f0.items()) for t in np.linspace(0.6,1,41)] for s in (0.05,0.15,0.3)],dtype=float)
                print(j,nm,q.degrees(),[ [round(0.6+0.01*i,2) for i in range(40) if np.sign(r[i])!=np.sign(r[i+1])] for r in G])

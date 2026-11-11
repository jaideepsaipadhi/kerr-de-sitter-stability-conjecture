import sys
sys.argv=['x','margin/sca_res.pkl','CQ2,CH2,CH1','scC','0','40','EA']
exec(open('gcover.py').read().split("todo=[tuple")[0])
for k,v in root.items():
    for q in v:
        ex=set((int(e[0]),int(e[2])) for e in q.to_dict())
        # lower-left Newton polygon vertices
        pts=sorted(ex); hull=[]
        for p in pts:
            if all(not(o[0]<=p[0] and o[1]<=p[1] and o!=p) for o in ex): hull.append(p)
        print(k,sorted(hull))
import numpy as np
for k,v in root.items():
    for q in v:
        d=q.to_dict(); lo=min(int(e[0])+int(e[2]) for e in d)
        h=[(int(e[0]),int(e[1]),int(e[2]),int(c)) for e,c in d.items() if int(e[0])+int(e[2])==lo]
        vals=[sum(c*g**a*b**bb*1**t for a,bb,t,c in h) for g in np.linspace(0,1,21) for b in np.linspace(0,1,21)]
        print(k,'lead deg',lo,'sign range',np.sign(min(vals)),np.sign(max(vals)), 'at g=0,b=0..', [np.sign(sum(c*0**a*b**bb for a,bb,t,c in h)) for b in (0,0.5,1)])

import sys
sys.argv=['x','margin/sca_res.pkl','CQ2,CH2,CH1','scC','0','40',sys.argv[1]]
src=open('gcover.py').read().split("import json")[0]+open('gcover.py').read().split("import json")[1].split("PK=f")[0]
exec(open('gcover.py').read().split("todo=[tuple")[0])
import numpy as np
B=np.linspace(0.01,0.99,30); Tt=np.linspace(0.01,0.99,30)
for k,v in root.items():
    out=[]
    for q in v:
        f0={(int(e[1]),int(e[2])):int(c) for e,c in q.to_dict().items() if int(e[0])==0}
        if not f0: out.append('ZERO'); continue
        vals=np.array([[sum(c*b**e[0]*t**e[1] for e,c in f0.items()) for t in Tt] for b in B])
        out.append((int(np.sign(vals.min())),int(np.sign(vals.max()))))
    print(k,S[k][0],out)
for k,v in root.items():
    for q in v:
        f0={(int(e[1]),int(e[2])):int(c) for e,c in q.to_dict().items() if int(e[0])==0}
        e00=lambda b,t: sum(c*b**e[0]*t**e[1] for e,c in f0.items())
        print(k, [float(np.sign(e00(b,t))) for b,t in [(0,0.5),(1,0.5),(0.5,0),(0.5,1),(0,0),(1,1),(0,1),(1,0)]])

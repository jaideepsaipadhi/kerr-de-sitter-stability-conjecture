import sys; sys.argv=['x','EA','0','0']
exec(open('cover2d.py').read().split("t0=time.time(); st=")[0])
import numpy as np
Bs=np.linspace(0,1,21)
for key,(s0,cs) in S.items():
    if key[0]!='H2': continue
    sg=[]
    for b in Bs:
        v=s0
        for q in cs:
            f0=sum(int(c)*b**int(e[1]) for e,c in q.to_dict().items() if int(e[0])==0); v*=np.sign(f0)
        sg.append(int(v))
    print(key, sg)

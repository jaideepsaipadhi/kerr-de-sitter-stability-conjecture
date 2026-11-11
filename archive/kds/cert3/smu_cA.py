exec(open('smu_corner.py').read().split("ch=sys.argv[1]")[0])
import numpy as np
C={k:[chart(q,'A')[0] for q in v[1]] for k,v in S.items()}
for k,v in C.items():
    out=[]
    for q in v:
        d=q.to_dict(); f0={(int(e[1]),int(e[2])):int(c) for e,c in d.items() if int(e[0])==0}
        B=np.linspace(0,1,41); G=np.linspace(0,1,41)
        vals=np.array([[sum(c*b**e[0]*g**e[1] for e,c in f0.items()) for g in G] for b in B])
        out.append((np.sign(vals.min()),np.sign(vals.max()), 'zero@' + str([(B[i],G[j]) for i,j in np.argwhere(np.abs(vals)<1e-9*np.abs(vals).max())][:3])))
    print(k,out)
def f0(q):
    d=q.to_dict(); return {(int(e[1]),int(e[2])):int(c) for e,c in d.items() if int(e[0])==0}
ev=lambda f,b,g: sum(c*b**e[0]*g**e[1] for e,c in f.items())
def prodsign(k,nm,b,g):
    s=S[(k,nm)][0]
    for q in C[(k,nm)]: s*=np.sign(ev(f0(q),b,g))
    return s
B=np.linspace(0.0125,0.9875,40)
neg=[(b,g) for b in B for g in B if prodsign('H2','P',b,g)<0]
print('H2-P<0 region pts',len(neg), 'beta range',min(x[0] for x in neg) if neg else None, max(x[0] for x in neg) if neg else None,'gamma range',min(x[1] for x in neg) if neg else None,max(x[1] for x in neg) if neg else None)
from collections import Counter
print(Counter((prodsign('G','F1',b,g),prodsign('G','Y',b,g),prodsign('G','P',b,g)) for b,g in neg))

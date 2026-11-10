import sys
sys.argv=['x','margin/sca_res.pkl','CQ2,CH2,CH1','scC','0','40','EA2']
exec(open('gcover.py').read().split("todo=[tuple")[0])
import numpy as np
D=48.0
def lead(cs,s0):
    hs=[]
    for q in cs:
        d=q.to_dict(); hs.append({(int(e[1]),int(e[2])):int(c) for e,c in d.items() if int(e[0])==0})
    def f(b,g):
        v=float(s0)
        for h in hs: v*=sum(c*b**e[0]*g**e[1] for e,c in h.items())
        return v
    return f
G=np.linspace(0,1,41); Bb=np.linspace(0,1,21)
for k in ('CQ2','CH2','CH1'):
    F1=lead(S[(k,'F1')][1] if False else root[(k,'F1')],S[(k,'F1')][0]); Y=lead(root[(k,'Y')],S[(k,'Y')][0]); P=lead(root[(k,'P')],S[(k,'P')][0])
    vals=np.array([[ (F1(b,g)*np.sqrt(D)-Y(b,g))/2 for g in G] for b in Bb])
    print(k,'lead value sign range',np.sign(vals.min()),np.sign(vals.max()),'min at',np.unravel_index(vals.argmin(),vals.shape), 'scale',np.abs(vals).max())

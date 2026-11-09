import numpy as np, pickle, sys, flint
exec(open('win_bern.py').read().split("which=")[0])
F1,Y,P=res['Q2']
c,fl=ctx.from_dict(F1).factor()
f=[f for f,m in fl if sum(f.degrees())>8][0]; A=to_arr(f)
def ev(x,v,t):
    i,j,k=np.indices(A.shape); return float(sum(int(A[a,b_,c_])*x**a*v**b_*t**c_ for a,b_,c_ in zip(i.flat,j.flat,k.flat) if A[a,b_,c_]!=0))
g=np.linspace(0,1,6)
for name,fn in [('r0=0',lambda p,q: ev(0,p,q)),('r0=1',lambda p,q: ev(1,p,q)),('v=0',lambda p,q: ev(p,0,q)),('v=1',lambda p,q: ev(p,1,q)),('tau=0',lambda p,q: ev(p,q,0)),('tau=1',lambda p,q: ev(p,q,1))]:
    vals=[[fn(p,q) for q in g] for p in g]
    print(name, np.sign(np.array(vals)).astype(int).tolist())

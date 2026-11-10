import sympy as sp, numpy as np, pickle
from fractions import Fraction as Fr
t,z0,A,be=sp.symbols('t z0 A beta'); T=pickle.load(open('T134.pkl','rb'))
P=[sp.Poly(Tk,t,z0,A,be) for Tk in T]
def ev(p,vals): return sum(Fr(int(c))*np.prod([v**int(e) for v,e in zip(vals,m)]) for m,c in zip(p.monoms(),p.coeffs()))
f=sp.lambdify((t,z0,A,be),T)
rng=np.random.default_rng(8); viol=0; real=0
for _ in range(200000):
    tt=1+10**rng.uniform(-2,3); zz=1+(tt-1)*rng.uniform(1e-2,1-1e-2); x0=zz-1
    a0=max(0.0,zz*((tt-1)/(2*x0**2)-3))
    AA=a0+10**rng.uniform(-3,4)*rng.integers(0,2)*(a0+1); need=(tt-1)/(2*x0)-3*x0-AA*x0/zz
    bb=max(need,0)+10**rng.uniform(-3,3)*(1+abs(need))
    c=np.array(f(tt,zz,AA,bb),float); sc=np.abs(c).max()
    if any(c[j]<-1e-9*sc and c[j+1]>1e-9*sc for j in range(1,6)):
        viol+=1; vals=[Fr(x) for x in (tt,zz,AA,bb)]
        ce=[ev(p,vals) for p in P]; sg=[(x>0)-(x<0) for x in ce if x!=0]
        n=sum(1 for i in range(len(sg)-1) if sg[i]!=sg[i+1])
        if n>1: real+=1; print([float(x) for x in vals],sg)
print('float viol',viol,'exact multi-change',real)

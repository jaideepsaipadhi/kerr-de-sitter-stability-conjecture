import sympy as sp, numpy as np, pickle
from fractions import Fraction as Fr
t,z0=sp.symbols('t z0'); G=pickle.load(open('G135.pkl','rb'))
rng=np.random.default_rng(1); 
for key in [('I',2),('I',3),('I',4)]:
    a,b,g1,g2=G[key][:4]; F=sp.lambdify((t,z0),[a,b,g1,g2]); real=0; fl=0
    for _ in range(60000):
        tt=1+10**rng.uniform(-2,3); zz=1+(tt-1)*rng.uniform(1e-3,1-1e-3); x=zz-1
        if (tt-1)/(2*x)-3*x<0: continue
        va,vb,v1,v2=F(tt,zz)
        if vb>=0 and (v1<=0 or v2<=0):
            fl+=1; s={t:sp.Rational(Fr(tt)),z0:sp.Rational(Fr(zz))}
            eb,e1,e2=[e.subs(s) for e in (b,g1,g2)]
            if eb>=0 and (e1<=0 or e2<=0): real+=1; print(key,tt,zz,float(eb),float(e1),float(e2))
    print(key,'float',fl,'real',real)

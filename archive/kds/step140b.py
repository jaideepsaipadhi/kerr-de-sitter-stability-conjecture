import sympy as sp, numpy as np, pickle, sys
t,z0,A,S,be,sg=sp.symbols('t z0 A S beta sigma',positive=True)
beta_star,cb,N0,Gn,Gd=pickle.load(open('P140.pkl','rb'))
fb=sp.lambdify((t,z0,A,S),beta_star); fG=sp.lambdify((t,z0,A,S,be,sg),Gn)
CC,DD=float(sys.argv[1]),float(sys.argv[2])
rng=np.random.default_rng(6); tot=bad=0; ex=[]
sgv=np.linspace(1e-4,1,300)
for i in range(400000):
    tt=1+10**rng.uniform(-2,3); zz=1+(tt-1)*rng.uniform(1e-2,1-1e-2); x0=zz-1
    a0=max(0.0,zz*((tt-1)/(2*x0**2)-3))
    AA=a0+10**rng.uniform(-3,4)*rng.integers(0,2)*(a0+1)
    SS=x0*rng.uniform(1e-3,1-1e-3); b=fb(tt,zz,AA,SS)
    if not (b>=CC*AA/zz+DD*zz): continue
    BB=3*zz-AA/zz+b
    if 2*x0*(AA+BB-3)<=tt-1: continue
    tot+=1; g=fG(tt,zz,AA,SS,b,sgv)
    if g.max()>0: bad+=1; ex.append((tt,x0/(tt-1),AA,SS/x0,b))
print(CC,DD,tot,bad,ex[:3])

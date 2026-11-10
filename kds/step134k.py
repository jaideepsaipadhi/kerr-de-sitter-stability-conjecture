import sympy as sp, numpy as np, pickle
t,z0,A,be=sp.symbols('t z0 A beta'); T=pickle.load(open('T134.pkl','rb'))
f=sp.lambdify((t,z0,A,be),T)
rng=np.random.default_rng(7); cnt={}; ex=[]
for _ in range(100000):
    tt=1+10**rng.uniform(-2,3); zz=1+(tt-1)*rng.uniform(1e-2,1-1e-2); x0=zz-1
    a0=max(0.0,zz*((tt-1)/(2*x0**2)-3))
    # A>=a0, beta>=0, p>1: A x0/z0+beta > (t-1)/(2x0)-3x0
    AA=a0+10**rng.uniform(-4,4)*rng.integers(0,2)*(a0+1); need=(tt-1)/(2*x0)-3*x0-AA*x0/zz
    bb=max(need,0)+10**rng.uniform(-5,3)
    c=np.array(f(tt,zz,AA,bb),float); sg=np.sign(c[c!=0]); ch=int(np.sum(sg[1:]!=sg[:-1]))
    cnt[ch]=cnt.get(ch,0)+1
    if ch!=1 and len(ex)<4: ex.append((tt,zz,AA,bb,np.sign(c)))
print(cnt,ex)

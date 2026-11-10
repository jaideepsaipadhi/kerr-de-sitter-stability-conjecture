import sympy as sp, numpy as np, pickle
t,z0,A,be,v=sp.symbols('t z0 A beta v'); T=pickle.load(open('T134.pkl','rb'))
f=[sp.lambdify((t,z0,A,be),Tk) for Tk in T]
rng=np.random.default_rng(5); shown=0
for _ in range(2000):
    tt=1+10**rng.uniform(-2,3); zz=1+(tt-1)*rng.uniform(1e-2,1-1e-2); x0=zz-1
    a0=max(0.0,zz*((tt-1)/(2*x0**2)-3))
    c=[f[k](tt,zz,a0,0.0) for k in range(7)]
    if any(c[k]>0 for k in (4,5,6)):
        r=np.roots(c[::-1]); pos=[x for x in r if abs(x.imag)<1e-9 and x.real>0]
        print(round(tt,3),round(zz,3),round(a0,2),np.sign(c).astype(int),len(pos)); shown+=1
        if shown>8: break

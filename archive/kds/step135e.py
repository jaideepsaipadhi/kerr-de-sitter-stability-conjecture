import sympy as sp, numpy as np, pickle
t,z0=sp.symbols('t z0'); G=pickle.load(open('G135.pkl','rb'))
rng=np.random.default_rng(2)
for key,val in G.items():
    a,b,g1,g2=val[:4]; F=sp.lambdify((t,z0),[g1,g2]); neg=[0,0]; n=0
    for _ in range(30000):
        tt=1+10**rng.uniform(-2,3); zz=1+(tt-1)*rng.uniform(1e-3,1-1e-3); x=zz-1
        if (key[0]=='I')!=((tt-1)/(2*x)-3*x>=0): continue
        n+=1; v=F(tt,zz); neg[0]+=v[0]<=0; neg[1]+=v[1]<=0
    print(key,n,neg)

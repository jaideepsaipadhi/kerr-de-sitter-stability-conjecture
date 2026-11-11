import sympy as sp, numpy as np, pickle, sys
t,z0,A,S,be,sg=sp.symbols('t z0 A S beta sigma',positive=True)
beta_star,cb,N0,Gn,Gd=pickle.load(open('P140.pkl','rb'))
m=A/(2*z0)+z0
g=sp.expand(N0+m*cb); Gm=sp.expand(Gn.subs(be,m))
pickle.dump((g,Gm),open('P140c.pkl','wb'))
fg=sp.lambdify((t,z0,A,S),g); fG=sp.lambdify((t,z0,A,S,sg),Gm)
rng=np.random.default_rng(7); tot=bad=0; badR=0; sgv=np.linspace(1e-4,1,300)
for i in range(40000):
    tt=1+10**rng.uniform(-2,3); zz=1+(tt-1)*rng.uniform(1e-3,1-1e-3); x0=zz-1
    AA=10**rng.uniform(-3,6); SS=x0*rng.uniform(1e-3,1-1e-3)
    if fg(tt,zz,AA,SS)<0: continue
    tot+=1; gg=fG(tt,zz,AA,SS,sgv).max()
    if gg>0:
        bad+=1
        a0=max(0.0,zz*((tt-1)/(2*x0**2)-3))
        if AA>=a0: badR+=1
print(tot,'bad (no R3 gate)',bad,'bad with A>=a0',badR)

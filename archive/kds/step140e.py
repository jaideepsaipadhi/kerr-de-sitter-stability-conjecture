import sympy as sp, numpy as np, pickle
t,z0,A,S,be,sg=sp.symbols('t z0 A S beta sigma',positive=True); ga=sp.Symbol('gamma')
beta_star,cb,N0,Gn,Gd=pickle.load(open('P140.pkl','rb'))
bsub=A/(2*z0)+z0+ga
Nfull=sp.expand(N0+cb*bsub)          # N(S) with beta = m + gamma
cA=Nfull.coeff(A); r=sp.expand(Nfull-cA*A)
Astar=sp.factor(-r/cA)
print('cA',sp.factor(cA))
Gsub=sp.together(Gn.subs(be,bsub).subs(A,Astar))
num,den=sp.fraction(Gsub)
print('den',sp.factor(den))
pickle.dump((Astar,num,den),open('P140e.pkl','wb'))
fA=sp.lambdify((t,z0,S,ga),Astar); fN=sp.lambdify((t,z0,S,ga,sg),num); fD=sp.lambdify((t,z0,S,ga,sg),den)
rng=np.random.default_rng(8); tot=badpos=badneg=nneg=0; sgv=np.linspace(1e-4,1,300)
for i in range(30000):
    tt=1+10**rng.uniform(-2,3); zz=1+(tt-1)*rng.uniform(1e-3,1-1e-3); x0=zz-1
    SS=x0*rng.uniform(1e-3,1-1e-3); gg=10**rng.uniform(-4,4)
    a=fA(tt,zz,SS,gg); v=(fN(tt,zz,SS,gg,sgv)/fD(tt,zz,SS,gg,sgv)).max()
    if a>=0: tot+=1; badpos+=v>0
    else: nneg+=1; badneg+=v>0
print('A>=0:',tot,'bad',badpos,' A<0:',nneg,'bad',badneg)
rng=np.random.default_rng(8)
for i in range(30000):
    tt=1+10**rng.uniform(-2,3); zz=1+(tt-1)*rng.uniform(1e-3,1-1e-3); x0=zz-1
    SS=x0*rng.uniform(1e-3,1-1e-3); gg=10**rng.uniform(-4,4)
    a=fA(tt,zz,SS,gg); vv=fN(tt,zz,SS,gg,sgv)/fD(tt,zz,SS,gg,sgv)
    if a>=0 and vv.max()>0: print(tt,(zz-1)/(tt-1),SS/x0,gg,a,vv.max(),np.abs(vv).max(),sgv[vv.argmax()])

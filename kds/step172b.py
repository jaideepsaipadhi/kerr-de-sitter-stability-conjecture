import sympy as sp, pickle, numpy as np
z,t,z0,S=sp.symbols('z t z0 S')
g0,g1,q0,q1,Qd=pickle.load(open('P172.pkl','rb'))
f=sp.lambdify((t,z0,S,z),[g0,g1,q0,q1])
rng=np.random.default_rng(1); n=2000000
x0=10**rng.uniform(-4,4,n); zz=1+x0; tt=zz+x0*10**rng.uniform(-6,7,n); Sv=x0*rng.uniform(1e-6,1-1e-6,n); zv=zz-Sv*rng.uniform(1e-9,1,n)
G0,G1,Q0,Q1=f(tt,zz,Sv,zv)
print('q1<0 frac',(Q1<0).mean(),' g1>0 frac',(G1>0).mean())
print('(iii) g1>0 & q1<0:',((G1>0)&(Q1<0)).sum())
print('(i) g0>=0 & q0<=0:',((G0>=0)&(Q0<=0)).sum())
As=-G0/G1; inr=(As>0)
QA=Q0+Q1*As
print('(ii) A* in range & Q(A*)<=0:',(inr&(QA<=0)).sum(), ' in range count',inr.sum())
W=Q1*G0-Q0*G1
print('ungated: q0<=0',(Q0<=0).sum(),' W<=0',(W<=0).sum(), ' W<=0 & g0>=0',((W<=0)&(G0>=0)).sum())

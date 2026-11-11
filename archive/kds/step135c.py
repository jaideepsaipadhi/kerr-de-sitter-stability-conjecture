import sympy as sp, numpy as np, pickle
t,z0,A,be=sp.symbols('t z0 A beta'); T=pickle.load(open('T134.pkl','rb'))
x0=z0-1; K=(t-1)/(2*x0)-3*x0; Ac2=z0*K/x0
G={}
for case,a0 in (('I',Ac2),('II',sp.Integer(0))):
  for j in (2,3,4):
    Tj,Tk=T[j],T[j+1]
    a=sp.together(Tj.subs({A:a0,be:0})); b=sp.together(Tk.subs({A:a0,be:0}))
    pA,pB=Tj.coeff(A),Tj.coeff(be); qA,qB=Tk.coeff(A),Tk.coeff(be)
    G[(case,j)]=(a,b,sp.together(-a*qA+pA*b),sp.together(-a*qB+pB*b),pA,pB,qA,qB)
pickle.dump(G,open('G135.pkl','wb'))
rng=np.random.default_rng(1)
for key,(a,b,g1,g2,pA,pB,qA,qB) in G.items():
    F=sp.lambdify((t,z0),[a,b,g1,g2,pA,pB,qA,qB]); bad=0; n=0
    for _ in range(60000):
        tt=1+10**rng.uniform(-2,3); zz=1+(tt-1)*rng.uniform(1e-3,1-1e-3); x=zz-1
        Kv=(tt-1)/(2*x)-3*x
        if (key[0]=='I')!=(Kv>=0): continue
        n+=1; va,vb,v1,v2,*cf=F(tt,zz)
        pass
        if vb>=0 and (v1<=0 or v2<=0): bad+=1
    print(key,n,'bad',bad)

import sympy as sp, pickle, flint
z,t,z0,S=sp.symbols('z t z0 S')
X,T,sg,rh=sp.symbols('X T sigma rho')
g0,g1,q0,q1,Qd=pickle.load(open('P172.pkl','rb'))
x0=X/(1-X); sub={z0:1+x0, t:1+x0+x0*T/(1-T), S:x0*sg}
def tocube(e,withz=True):
    e2=e.subs(z, z0-S*rh) if withz else e
    e2=sp.together(e2.subs(sub))
    n,d=sp.fraction(sp.factor(e2))
    return n,sp.factor(d)
out={}
for nm,e in (('g0',g0),('g1',g1),('q0',q0),('W',sp.expand(q1*g0-q0*g1))):
    n,d=tocube(e); out[nm]=(sp.Poly(sp.expand(n),X,T,sg,rh),d)
    print(nm,'deg',[out[nm][0].degree(v) for v in (X,T,sg,rh)],'den',d,flush=True)
pickle.dump(out,open('P173.pkl','wb'))

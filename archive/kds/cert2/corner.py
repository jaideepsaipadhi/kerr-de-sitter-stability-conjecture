import sympy as sp, pickle
exec(open('a0_bern.py').read().split("from math import comb")[0].replace("print(nm,'deg'","pass #"))
out=pickle.load(open('a0_scaled.pkl','rb')); nu,_=out['pos']
P=sp.expand(nu.subs(be,1/bt)*bt**sp.Poly(nu,T,ze,y,be).degree(be))
eps=sp.Symbol('eps',positive=True)
# scale T=a*eps, bt=b*eps, 1-y=c*eps
a,b,c=sp.symbols('a b c',positive=True)
Ps=sp.expand(P.subs({T:a*eps,bt:b*eps,y:1-c*eps}))
lead=sp.Poly(Ps,eps); low=min(m[0] for m in lead.monoms())
print('lowest eps order',low); L=sp.factor(lead.coeff_monomial(eps**low)); print('leading form',L)

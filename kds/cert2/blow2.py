import pickle, sympy as sp
T,ze,y,om=sp.symbols('T ze y om',positive=True); a,c,eps=sp.symbols('a c eps',positive=True)
res=pickle.load(open('a0w.pkl','rb'))
for nm in ('pos','neg'):
    P=sp.expand(res[nm][0].as_expr().subs({ze:1-a*eps,y:1-c*eps}))
    Pe=sp.Poly(P,eps); low=min(m[0] for m in Pe.monoms()); Lf=sp.factor(Pe.coeff_monomial(eps**low))
    print(nm,'order',low,'leading',Lf)
for nm in ('pos','neg'):
    P=sp.expand(res[nm][0].as_expr().subs({ze:1,y:1-c*eps}))
    Pe=sp.Poly(P,eps); low=min(m[0] for m in Pe.monoms()); print(nm,'on ze=1: order',low,sp.factor(Pe.coeff_monomial(eps**low)))

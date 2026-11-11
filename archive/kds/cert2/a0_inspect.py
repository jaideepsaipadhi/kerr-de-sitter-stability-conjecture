import sympy as sp, pickle
T,ze,y,be=sp.symbols('T ze y be',positive=True); bt=sp.Symbol('bt',positive=True)
out=pickle.load(open('a0_scaled.pkl','rb'))
for nm,(nu,dn) in out.items():
    P=sp.Poly(nu,T,ze,y,be); db=P.degree(be)
    Q=sp.Poly(sp.expand(nu.subs(be,1/bt)*bt**db),T,ze,y,bt)
    print(nm,'deg T,ze,y,bt',Q.degree(T),Q.degree(ze),Q.degree(y),Q.degree(bt))
    for var in (T,ze,y,bt):
        e=sp.factor(Q.as_expr().subs(var,0)); print('   at',var,'=0:',str(e)[:160])
    pickle.dump(Q,open(f'a0_{nm}.pkl','wb'))

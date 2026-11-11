import pickle, sympy as sp, sys
T,ze,y,p=sp.symbols('T ze y p',positive=True)
res=pickle.load(open(f'surf_{sys.argv[1]}.pkl','rb'))
for k in ('pos','neg'):
    P=res[k][0]; mz=min(m[1] for m in P.monoms()); print(k,'min power of ze',mz)
    Q=sp.Poly(sp.expand(P.as_expr()/ze**mz),T,ze,y,p)
    print('  at ze=0:',sp.factor(Q.as_expr().subs(ze,0)))

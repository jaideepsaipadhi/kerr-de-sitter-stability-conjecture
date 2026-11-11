import pickle, sympy as sp
T,ze,y,om=sp.symbols('T ze y om',positive=True)
res=pickle.load(open('a0w.pkl','rb'))
for nm in ('pos','neg','floor','pp','mup'):
    P=res[nm][0].as_expr()
    for (Tv,e1,e2) in [(sp.Rational(6,7),sp.Rational(1,1000),sp.Rational(1,1000)),(sp.Rational(6,7),sp.Rational(1,10**5),sp.Rational(1,10**5)),(sp.Rational(5,6),sp.Rational(1,1000),sp.Rational(1,1000)),(sp.Rational(9,10),sp.Rational(1,1000),sp.Rational(1,1000))]:
        print(nm,float(Tv),float(e1),sp.N(P.subs({T:Tv,ze:1-e1,y:1-e2,om:1}),6))

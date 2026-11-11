import pickle, sympy as sp
T,ze,c,e=sp.symbols('T ze c e',positive=True)
out=pickle.load(open('/home/claude/kds/p1_2d.pkl','rb'))
for k in ('Ac2','G0_I','al4_I','al5_I','al6_I','al3_I','al2_I','G0_II','al4_II','al3_II'):
    P,dn=out[k]
    for (a,b) in ((1,1),(1,2),(2,1)):  # T = e^a, ze = c e^b
        ex=sp.expand(P.as_expr().subs({T:e**a,ze:c*e**b}))
        pe=sp.Poly(ex,e); lo=min(m[0] for m in pe.monoms()); lead=sp.factor(pe.coeff_monomial(e**lo))
        dd=sp.factor(dn.subs({T:e**a,ze:c*e**b}))
        print(k,(a,b),'lead',lead,'| den',str(dd)[:40])

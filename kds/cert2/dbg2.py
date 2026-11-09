import sys; sys.argv=['p1chart.py','R2']
src=open('p1chart.py').read().split("for c in ('I','II'):")[0]
exec(src)
b=BB
for k in ['Ac2','G0_II','al4_II','al5_II','al6_II','al3_II','al2_II','v2_II','v3_II']:
    print(k,'min',float(mn(b[k])),'max',float(mx(b[k])))
import sympy as sp
P,dn=out['Ac2']; print(sp.N((P.as_expr()/dn).subs({T:sp.Rational(1,10),ze:sp.Rational(1,10)})))

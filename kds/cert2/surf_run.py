import pickle, sympy as sp, sys, random
from bern_core import *
T,ze,y,p=sp.symbols('T ze y p',positive=True); V=(T,ze,y,p)
name=sys.argv[1]; lo=[float(a) for a in sys.argv[2:6]]; hi=[float(a) for a in sys.argv[6:10]]
res=pickle.load(open(f'surf_{name}.pkl','rb'))
# sign of each denominator: evaluate at random interior points; require constant sign (factors of fixed sign)
sg={}
for k,(P,dn) in res.items():
    vals=set()
    for _ in range(20):
        pt={v:sp.Rational(random.randint(1,999),1000) for v in V}
        dv=dn.subs(pt) if hasattr(dn,'subs') else dn
        vals.add(sp.sign(dv))
    assert len(vals)==1, (k,vals); sg[k]=int(list(vals)[0])
print('den signs',sg)
sub={v:sp.Rational(str(lo[i]))+(sp.Rational(str(hi[i]))-sp.Rational(str(lo[i])))*v for i,v in enumerate(V)}
BB={}
for k,(P,dn) in res.items():
    if P.is_zero: continue
    Q=sp.Poly(sp.expand(P.as_expr().subs(sub,simultaneous=True)*sg[k]),*V)   # now: quantity has the sign of Q
    BB[k]=to_bern(Q,V)
def rule(b):
    # constraints Q >= 0 required (pp strictly); violated if max < 0
    for k in ('bulk','c2','floor','Apos'):
        if k in b and mx(b[k])<0: return True
    if 'pp' in b and mx(b['pp'])<=0: return True
    # F <= 0  <=>  Q_F <= 0
    return (mx(b['mup'])<=0 or mx(b['pos'])<=0) and (mn(b['mup'])>=0 or mx(b['neg'])<=0)
r,deep=bb(BB,rule,4,maxdep=44)
print(r); print([tuple((round(lo[i]+(hi[i]-lo[i])*a,5),round(lo[i]+(hi[i]-lo[i])*c,5)) for i,(a,c) in enumerate(d)) for d in deep[:4]])

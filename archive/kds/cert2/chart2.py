import pickle, sympy as sp, sys
from bern_core import *
T,ze,y,be=sp.symbols('T ze y be',positive=True); bt=sp.Symbol('bt',positive=True); V=(T,ze,y,bt)
exec(open('a0_bern.py').read().split("from math import comb")[0].replace("print(nm,'deg'","pass #"))
out=pickle.load(open('a0_scaled.pkl','rb'))
polys={}
for nm in ('pos','neg'):
    nu,_=out[nm]; P=sp.Poly(nu,T,ze,y,be); db=P.degree(be); polys[nm]=sp.Poly(sp.expand(nu.subs(be,1/bt)*bt**db),*V)
for nm in ('c2','pp','bulk','mup'): polys[nm]=cons[nm][0]
lo=[float(a) for a in sys.argv[1:5]]; hi=[float(a) for a in sys.argv[5:9]]
# restrict polynomials to the sub-box by affine substitution
sub={v:lo[i]+(hi[i]-lo[i])*v for i,v in enumerate(V)}
sub={v:sp.Rational(str(lo[i]))+(sp.Rational(str(hi[i]))-sp.Rational(str(lo[i])))*v for i,v in enumerate(V)}
BB={nm:to_bern(sp.Poly(sp.expand(P.as_expr().subs(sub,simultaneous=True)),*V),V) for nm,P in polys.items()}
def rule(b):
    if mx(b['c2'])<0 or mx(b['pp'])<=0 or mx(b['bulk'])<0: return True
    return (mx(b['mup'])<=0 or mn(b['pos'])>=0) and (mn(b['mup'])>=0 or mn(b['neg'])>=0)
r,deep=bb(BB,rule,4,maxdep=48)
print(r); print([tuple((round(lo[i]+(hi[i]-lo[i])*a,5),round(lo[i]+(hi[i]-lo[i])*c,5)) for i,(a,c) in enumerate(d)) for d in deep[:4]])

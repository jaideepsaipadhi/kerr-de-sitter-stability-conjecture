import pickle, sympy as sp
from bern_core import *
T,ze,y,om=sp.symbols('T ze y om',positive=True); V=(T,ze,y,om)
res=pickle.load(open('a0w.pkl','rb'))
BB={nm:to_bern(P,V) for nm,(P,dn) in res.items()}
def rule(b):
    if mx(b['floor'])<0 or mx(b['pp'])<=0 or mx(b['c2'])<0: return True
    okp=mx(b['mup'])<=0 or mn(b['pos'])>=0
    okn=mn(b['mup'])>=0 or mn(b['neg'])>=0
    return okp and okn
r,deep=bb(BB,rule,4,maxdep=60)
print(r); print([tuple((round(a,5),round(c,5)) for a,c in d) for d in deep[:6]])

import pickle, sympy as sp, random, sys
from bern_core import to_bern, bb, mn, mx
T,ze=sp.symbols('T ze',positive=True); V=(T,ze)
out=pickle.load(open('/home/claude/kds/p1_2d.pkl','rb'))
sg={}
for k,(P,dn) in out.items():
    s={int(sp.sign(dn.subs({T:sp.Rational(random.randint(1,999),1000),ze:sp.Rational(random.randint(1,999),1000)}))) for _ in range(10)}
    assert len(s)==1; sg[k]=s.pop()
def strip(P):
    e=sp.factor(P.as_expr()); rem=sp.Integer(1); keep=sp.Integer(1)
    for f_,m_ in sp.factor_list(P.as_expr())[1]:
        if f_ in (T,ze,1-T,T-1,1-ze,ze-1):
            # positive on the open box up to sign: (T-1),(ze-1) negative
            if f_ in (T-1,ze-1): keep*=(-1)**m_
            continue
        rem*=f_**m_
    c=sp.factor_list(P.as_expr())[0]
    return sp.Poly(sp.expand(c*keep*rem),*V)
BB={}
for k,(P,dn) in out.items():
    Q=strip(sp.Poly(P.as_expr()*sg[k],*V)); BB[k]=to_bern(Q,V)   # each array has the SIGN of the quantity
def rule_case(c):
    def r(b):
        ac=b['Ac2']
        if c=='I' and mn(ac)>0: return True      # not in case I
        if c=='II' and mx(ac)<0: return True     # not in case II
        if not mn(b['G0_'+c])>=0: return False
        if not (mx(b['al4_'+c])<0 and mx(b['al5_'+c])<0 and mx(b['al6_'+c])<0): return False
        if mx(b['al3_'+c])<0: return True
        return mn(b['al2_'+c])>0 and mn(b['v2_'+c])>0 and mn(b['v3_'+c])>0
    return r
for c in ('I','II'):
    res,deep=bb(BB,rule_case(c),2,tmax=25,maxdep=40)
    print('case',c,res,[tuple((round(a,6),round(b,6)) for a,b in d) for d in deep[:4]])

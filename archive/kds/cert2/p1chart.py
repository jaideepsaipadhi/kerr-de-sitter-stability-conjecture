# (P1) 2-D certificate with an optional blow-up chart:  T = phi_T(u,v), ze = phi_z(u,v) on [0,1]^2 (scaled).
import pickle, sympy as sp, random, sys
from bern_core import to_bern, bb, mn, mx
T,ze,u,w=sp.symbols('T ze u w',positive=True)
out=pickle.load(open('/home/claude/kds/p1_2d.pkl','rb'))
chart=sys.argv[1]; d=sp.Rational(sys.argv[2]) if len(sys.argv)>2 else sp.Rational(1)
maps={'orig':(u,w),
      'R1':(sp.Rational(1,10)+sp.Rational(9,10)*u,w), 'R2':(u/10,sp.Rational(1,10)+sp.Rational(9,10)*w),
      'A':(d*u,d*u*w),        # zeta = T*w  (zeta <= T), T in [0,d]
      'B':(d*u*w,d*u),        # T = zeta*w  (T <= zeta), zeta in [0,d]
      'A2':(d*u,d*u**2*w),    # zeta <= T^2
      'B2':(d*u**2*w,d*u)}    # T <= zeta^2
mT,mz=maps[chart]
def strip(e):
    e=sp.factor(e); c,fl=sp.factor_list(e); rem=sp.Integer(1); sgn=1
    for f_,m_ in fl:
        if f_ in (u,w): continue
        if f_ in (u-1,w-1): sgn*=(-1)**m_; continue
        rem*=f_**m_
    return sp.Poly(sp.expand(c*sgn*rem),u,w)
BB={}
for k,(P,dn) in out.items():
    ex=P.as_expr().subs({T:mT,ze:mz},simultaneous=True); de=dn.subs({T:mT,ze:mz},simultaneous=True)
    s={int(sp.sign(de.subs({u:sp.Rational(random.randint(1,999),1000),w:sp.Rational(random.randint(1,999),1000)}))) for _ in range(8)}
    assert len(s)==1,k; BB[k]=to_bern(strip(ex*s.pop()),(u,w))
def rule_case(c):
    def r(b):
        ac=b['Ac2']
        if c=='I' and mn(ac)>0: return True
        if c=='II' and mx(ac)<0: return True
        if not mn(b['G0_'+c])>=0: return False
        if not (mx(b['al4_'+c])<0 and mx(b['al5_'+c])<0 and mx(b['al6_'+c])<0): return False
        if mx(b['al3_'+c])<0: return True
        return mn(b['al2_'+c])>0 and mn(b['v2_'+c])>0 and mn(b['v3_'+c])>0
    return r
for c in ('I','II'):
    res,deep=bb(BB,rule_case(c),2,tmax=12,maxdep=36)
    print(chart,str(d),'case',c,res,[tuple((round(a,6),round(b,6)) for a,b in dd) for dd in deep[:3]])

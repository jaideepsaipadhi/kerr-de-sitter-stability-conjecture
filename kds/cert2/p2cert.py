import pickle, sympy as sp, random, sys
from bern_core import to_bern, bb, mn, mx
T,ze,sg=sp.symbols('T ze sg',positive=True); u,w,v=sp.symbols('u w v',positive=True)
out=pickle.load(open('/home/claude/kds/p2_3d.pkl','rb'))
chart=sys.argv[1]; d=sp.Rational(sys.argv[2]) if len(sys.argv)>2 else 1
dl=sp.Rational(1,10)
Tr={'R1':sp.Rational(1,10)+sp.Rational(9,10)*u,'R2':u/10}
def zmaps(Tm,zlo=0):
    # zeta in [zlo,1]; corner (zeta,sg)->(1,0) blown up
    span=1-zlo
    return {'m':(Tm,zlo+span*(1-dl)*w,v),'n':(Tm,1-dl*span+dl*span*w,dl+(1-dl)*v),
            'p':(Tm,1-dl*span*w,dl*span*w*v),'q':(Tm,1-dl*span*w*v,dl*span*w)}
maps={'orig':(u,w,v),'R1':(sp.Rational(1,10)+sp.Rational(9,10)*u,w,v),'R2':(u/10,sp.Rational(1,10)+sp.Rational(9,10)*w,v),
      'A':(d*u,d*u*w,v),'B':(d*u*w,d*u,v)}
for base,zlo in (('R1',0),('R2',sp.Rational(1,10))):
    for k_,m_ in zmaps(Tr[base],zlo).items(): maps[base+k_]=m_
# corner blow-up in (X, sigma) at (0, 1):  pieces for a chart's (x-variable, v)
def corner_pieces(xv):
    e=dl
    return {'_x':{xv:e+(1-e)*xv},                       # x >= e
            '_y':{xv:e*xv, v:(1-e)*v},                  # x <= e, sigma <= 1-e
            '_p':{xv:e*xv, v:1-e*xv*v},                 # x <= e, 1-sigma <= x
            '_q':{xv:e*xv*v, v:1-e*xv}}                 # 1-sigma <= e, x <= 1-sigma
base_names=list(maps.keys())
for bn in base_names:
    xv=w if bn=='B' else u
    for suf,subs in corner_pieces(xv).items():
        bm=maps[bn]; maps[bn+suf]=tuple(sp.expand(c.subs(subs,simultaneous=True)) for c in bm)
mT,mz,ms=maps[chart]
def strip(e):
    c,fl=sp.factor_list(sp.factor(e)); rem=sp.Integer(1); sgn=1
    for f_,m_ in fl:
        if f_ in (u,w,v): continue
        if f_ in (u-1,w-1,v-1): sgn*=(-1)**m_; continue
        rem*=f_**m_
    return sp.Poly(sp.expand(c*sgn*rem),u,w,v)
BB={}
for k,(P,dn) in out.items():
    ex=P.as_expr().subs({T:mT,ze:mz,sg:ms},simultaneous=True); de=dn.subs({T:mT,ze:mz,sg:ms},simultaneous=True)
    s={int(sp.sign(de.subs({u:sp.Rational(random.randint(1,999),1000),w:sp.Rational(random.randint(1,999),1000),v:sp.Rational(random.randint(1,999),1000)}))) for _ in range(8)}
    assert len(s)==1,k; BB[k]=to_bern(strip(ex*s.pop()),(u,w,v))
def rule(b):
    if not mx(b['BC'])<0: return False
    if mn(b['AK'])>=0: return mx(b['NAK'])<0
    if mx(b['AK'])<=0: return mx(b['NA0'])<0
    return mx(b['NAK'])<0 and mx(b['NA0'])<0
res,deep=bb(BB,rule,3,tmax=45,maxdep=40)
print(chart,str(d),res,[tuple((round(float(a),5),round(float(c),5)) for a,c in dd) for dd in deep[:3]])

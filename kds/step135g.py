import sympy as sp, pickle, sys, itertools
from nnd import cert
t,z0=sp.symbols('t z0'); u,R=sp.symbols('u R'); G=pickle.load(open('G135.pkl','rb')); x=sp.Symbol('x')
sub={'I': {t:1+6*x**2/u}, 'II': {t:1+x+(6*x**2-x)*u}}
xs={'I': u/6+R/(1-R), 'II': sp.Rational(1,6)+R/(1-R)}
def sign_factor(f):
    P=sp.Poly(f,u,R)
    if P.total_degree()<=1:
        vals=[f.subs({u:a,R:b}) for a,b in itertools.product((0,1),repeat=2)]
        if all(v>0 for v in vals) or (all(v>=0 for v in vals) and f.subs({u:sp.Rational(1,2),R:sp.Rational(1,2)})>0): return 1
        if all(v<0 for v in vals) or (all(v<=0 for v in vals) and f.subs({u:sp.Rational(1,2),R:sp.Rational(1,2)})<0): return -1
        return None
    if cert(f,[u,R],tmax=15)['ok']: return 1
    if cert(-f,[u,R],tmax=15)['ok']: return -1
    return None
def sign_expr(e,case):
    e=sp.together(e.subs(z0,1+x).subs(sub[case]).subs(x,xs[case]))
    n,d=sp.fraction(e); s=1
    for part in (n,d):
        c,fl=sp.factor_list(part); s*= 1 if c>0 else -1
        for f,m in fl:
            sg=sign_factor(f)
            if sg is None: return None
            if m%2: s*=sg
    return s
for case in sys.argv[1:]:
  for j in (2,3,4):
    a,b,g1,g2=G[(case,j)][:4]
    print(case,j,'g1',sign_expr(g1,case),'g2',sign_expr(g2,case),flush=True)

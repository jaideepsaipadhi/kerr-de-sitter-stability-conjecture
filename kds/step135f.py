import sympy as sp, pickle, sys
from nnd import cert
t,z0=sp.symbols('t z0'); u,R=sp.symbols('u R'); G=pickle.load(open('G135.pkl','rb'))
x=sp.Symbol('x')
sub={'I': {t:1+6*x**2/u}, 'II': {t:1+x+(6*x**2-x)*u}}
xs={'I': u/6+R/(1-R), 'II': sp.Rational(1,6)+R/(1-R)}
def signed(e,case):
    e=sp.together(e.subs(z0,1+x).subs(sub[case]).subs(x,xs[case]))
    n,d=sp.fraction(sp.factor(e)); return n,d
key=(sys.argv[1],int(sys.argv[2])); which=sys.argv[3]
a,b,g1,g2=G[key][:4]; e={'g1':g1,'g2':g2,'b':b}[which]
n,d=signed(e,key[0])
print(key,which,'num factors:',[ (sp.Poly(f,u,R).total_degree(),m) for f,m in sp.factor_list(n)[1]],'const',sp.factor_list(n)[0])
print('den:',sp.factor(d))
for f,m in sp.factor_list(n)[1]:
    if sp.Poly(f,u,R).total_degree()<=1: print(' small',f); continue
    r=cert(f,[u,R],tmax=20); r2=None if r['ok'] else cert(-f,[u,R],tmax=20)
    print(' factor deg',sp.Poly(f,u,R).degree(u),sp.Poly(f,u,R).degree(R),'mult',m,'pos' if r['ok'] else ('neg' if r2 and r2['ok'] else ('FAIL',r,r2)))

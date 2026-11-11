import pickle, sympy as sp, random
t_,z0_,A_,B_=sp.symbols('t z0 A B',positive=True)
Tj=pickle.load(open('/home/claude/kds/T118.pkl','rb'))
lam=2*(t_-z0_)/(7*t_-6*z0_-1); G=sp.together(Tj[2]-lam*Tj[3])
Gn,Gd=sp.fraction(G); GnP=sp.Poly(sp.expand(Gn),A_,B_)
print('G numerator A-coef:',sp.factor(GnP.coeff_monomial(A_)),' denominator:',sp.factor(Gd))
print('G numerator B-coef:',sp.factor(GnP.coeff_monomial(B_)))
fl=sp.Rational(11,2)*t_-sp.Rational(5,4)
targets={'T4(0,Bf)':Tj[4].subs({A_:0,B_:fl}),'T5(0,Bf)':Tj[5].subs({A_:0,B_:fl}),'T6(0,Bf)':Tj[6].subs({A_:0,B_:fl}),
         'Gnum(Bf)':GnP.as_expr().subs({A_:0,B_:fl})}
from fractions import Fraction
fs={k:sp.lambdify((t_,z0_),v,'sympy') for k,v in targets.items()}
mx={k:None for k in fs}
random.seed(1)
for _ in range(3000):
    tv=1+Fraction(10)**random.choice(range(-2,5))*Fraction(random.randint(1,1000),1000)
    z0v=1+Fraction(random.randint(1,999),1000)*(tv-1)
    for k,f in fs.items():
        val=f(sp.Rational(tv),sp.Rational(z0v)); s=sp.sign(val)
        mx[k]=(mx[k] or set())|{int(s)}
print('signs seen:',mx)

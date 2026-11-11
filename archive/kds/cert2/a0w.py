# A = 0 surface in variables (T, ze, y, om) with om = 10/w in (0,1]: the bulk constraint is a box side.
import sympy as sp, pickle, time
T,ze,y,om=sp.symbols('T ze y om',positive=True); V=(T,ze,y,om)
t=1/T; x0=ze*(t-1); z0=1+x0; x=y*x0; z=1+x; s=(1-y)*x0; mu0=z0*x0*(t-z0)
B=3*z0+10*mu0/(8*s**2*om)
L=B-3*z; g=3/L; mu=z*x*(t-z); mup=-3*z**2+2*(1+t)*z-t
v=1-y; q=sp.Rational(9,20)*v/(x*z); qp=sp.Rational(9,20)*(-(1/x0)/(x*z)-v*(2*z-1)/(x*z)**2)
res={}
for nm,c in (('pos',sp.Rational(17,10)),('neg',sp.Rational(17,11))):
    h=c*g+q; hu=sp.Rational(17,10)*g+q
    F=sp.Rational(1,2)*mup*h+sp.Rational(85,100)*mu*g**2+sp.Rational(1,2)*mu*hu**2+sp.Rational(1,2)*mu*qp-sp.Rational(5455,10000)*s*L*g-s*L*q
    nu,dn=sp.fraction(sp.factor(sp.together(F*T)))
    res[nm]=(sp.Poly(sp.expand(nu),*V),dn); print(nm,'deg',[res[nm][0].degree(v_) for v_ in V],'den',dn,flush=True)
fl=(B-(sp.Rational(11,2)*t-sp.Rational(5,4)))
nu,dn=sp.fraction(sp.factor(sp.together(fl))); res['floor']=(sp.Poly(sp.expand(nu),*V),dn); print('floor den',dn)
pp=2*x0*(B-3)-(t-1); nu,dn=sp.fraction(sp.factor(sp.together(pp))); res['pp']=(sp.Poly(sp.expand(nu),*V),dn); print('pp den',dn)
c2=6*ze**2*(1-T)-T; res['c2']=(sp.Poly(c2,*V),1)
nu,dn=sp.fraction(sp.factor(sp.together(mup))); res['mup']=(sp.Poly(sp.expand(nu),*V),dn); print('mup den',dn)
pickle.dump(res,open('a0w.pkl','wb'))

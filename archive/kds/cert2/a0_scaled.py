# B3 on the A = 0 surface, divided by t, as an explicit rational function of (T = 1/t, ze, y, be) -- for t -> infinity.
import sympy as sp, time
t0=time.time()
T,ze,y,be=sp.symbols('T ze y be',positive=True)
t=1/T
B=be*(sp.Rational(11,2)*t-sp.Rational(5,4))
x0=ze*(t-1); z0=1+x0; x=y*x0; z=1+x; s=(1-y)*x0
L=B-3*z; g=3/L
mu=z*x*(t-z); mup=-3*z**2+2*(1+t)*z-t
v=1-y; q=sp.Rational(9,20)*v/(x*z); qp=sp.Rational(9,20)*(-(1/x0)/(x*z)-v*(2*z-1)/(x*z)**2)
out={}
for nm,c in (('pos',sp.Rational(17,10)),('neg',sp.Rational(17,11))):
    h=c*g+q; hu=sp.Rational(17,10)*g+q
    F=sp.Rational(1,2)*mup*h+sp.Rational(85,100)*mu*g**2+sp.Rational(1,2)*mu*hu**2+sp.Rational(1,2)*mu*qp-sp.Rational(5455,10000)*s*L*g-s*L*q
    Ft=sp.factor(sp.together(F*T))
    nu,dn=sp.fraction(Ft)
    out[nm]=(sp.expand(nu),dn)
    print(nm,'num terms',len(sp.Add.make_args(sp.expand(nu))),'denom',dn, round(time.time()-t0,1),flush=True)
import pickle; pickle.dump(out,open('a0_scaled.pkl','wb'))

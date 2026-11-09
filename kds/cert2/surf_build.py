# Build B3 numerator polynomials and constraint polynomials on a parametrised surface (A(p), B(p)).
import sympy as sp, pickle, sys, time
T,ze,y,p=sp.symbols('T ze y p',positive=True); V=(T,ze,y,p)
def build(A,B,name):
    t0=time.time()
    t=1/T; x0=ze*(t-1); z0=1+x0; x=y*x0; z=1+x; s=(1-y)*x0; mu0=z0*x0*(t-z0)
    L=A/z+B-3*z; g=(A/z**2+3)/L; mu=z*x*(t-z); mup=-3*z**2+2*(1+t)*z-t
    v=1-y; q=sp.Rational(9,20)*v/(x*z); qp=sp.Rational(9,20)*(-(1/x0)/(x*z)-v*(2*z-1)/(x*z)**2)
    res={}
    for nm,c in (('pos',sp.Rational(17,10)),('neg',sp.Rational(17,11))):
        h=c*g+q; hu=sp.Rational(17,10)*g+q
        F=sp.Rational(1,2)*mup*h+sp.Rational(85,100)*mu*g**2+sp.Rational(1,2)*mu*hu**2+sp.Rational(1,2)*mu*qp-sp.Rational(5455,10000)*s*L*g-s*L*q
        nu,dn=sp.fraction(sp.factor(sp.together(F)))
        res[nm]=(sp.Poly(sp.expand(nu),*V),sp.factor(dn))
    L0=A/z0+B-3*z0
    for nm,ex in (('floor',B-(sp.Rational(11,2)*t-sp.Rational(5,4))),('bulk',8*L0*s**2-10*mu0),
                  ('pp',2*x0*(A+B-3)-(t-1)),('c2',2*x0**2*(A/z0+3)-(t-1)),('Apos',A),('mup',mup)):
        nu,dn=sp.fraction(sp.factor(sp.together(ex))); res[nm]=(sp.Poly(sp.expand(nu),*V),sp.factor(dn))
    for k,(P,dn) in res.items(): print(name,k,'deg',[P.degree(v_) for v_ in V],'den',str(dn)[:110],flush=True)
    pickle.dump(res,open(f'surf_{name}.pkl','wb')); print('built',round(time.time()-t0,1))
if __name__=='__main__':
    which=sys.argv[1]
    t=1/T; floor=sp.Rational(11,2)*t-sp.Rational(5,4)
    if which=='floor':   # B = floor, A = floor * p/(1-p), p in [0,1)
        build(floor*p/(1-p),floor,'floor')
    if which=='floorw':  # B = floor, A from w = 10/p  (p = omega in (0,1])
        x0=ze*(t-1); z0=1+x0; s=(1-y)*x0; mu0=z0*x0*(t-z0)
        build(z0*(10*mu0/(8*s**2*p)-floor+3*z0),floor,'floorw')

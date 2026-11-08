# Item 105: (C), (D) for the rational R3 multiplier as polynomial inequalities.
import sympy as sp, time
t0=time.time()
z,t,z0,A,B=sp.symbols('z t z0 A B',positive=True)
de,be,p=sp.Rational(7,10),8,sp.Rational(9,20)
mu=z*(z-1)*(t-z); mup=sp.diff(mu,z)
Lz=A/z+B-3*z; Lp=sp.diff(Lz,z)
g=-Lp/Lz
kap=(Lz/mu).subs(z,z0)
s=z0-z; x0=z0-1; u=(z-1)/x0
w=be*kap*s**2
h=(1+de)*g*w/(1+w)+p*(1-u)/((z-1)*z)
hp=sp.diff(h,z)
H=2*mu+s*mup
C=s*h*mu-H
D=sp.Rational(1,2)*(mup*h+mu*(hp+h**2))-s*(Lz*h+Lp)
for name,E in (('C',C),('D',D)):
    nu,dn=sp.fraction(sp.factor(sp.together(E)))
    P=sp.Poly(sp.expand(nu),z,t,z0,A,B)
    print(name,'num total deg',P.total_degree(),'terms',len(P.terms()),'| denom',sp.factor(dn), round(time.time()-t0,1),flush=True)

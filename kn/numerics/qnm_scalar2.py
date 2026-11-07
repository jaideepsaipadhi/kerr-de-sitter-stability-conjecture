import sympy as sp, numpy as np
from scipy.linalg import eig
u,w=sp.symbols('u w'); s,a,Q,m,lam,M=sp.symbols('s a Q m lam M')
def build(sv,av,Qv,mv,lamv,Mv=1):
    rp=Mv+sp.sqrt(Mv**2-av**2-Qv**2); rm=Mv-sp.sqrt(Mv**2-av**2-Qv**2)
    r=rp/(1-u); Dl=(r-rp)*(r-rm); K=(r**2+av**2)*w-av*mv
    OH=av/(rp**2+av**2); sig=(w-mv*OH)*(rp**2+av**2)/(rp-rm)
    A=-sv-sp.I*sig; B=-1-sv+sp.I*sig+2*sp.I*Mv*w
    h=sp.Function('h')(u)
    P=(r-rp)**A*(r-rm)**B*sp.exp(sp.I*w*r)
    R=P*h; dr=lambda f: sp.diff(f,u)*(1-u)**2/rp
    eq=Dl**2*dr(dr(R))+(sv+1)*Dl*sp.diff(Dl,u)*(1-u)**2/rp*dr(R)+(K**2-2*sp.I*sv*(r-Mv)*K+Dl*(4*sp.I*sv*w*r-lamv))*R
    eq=sp.expand(sp.simplify(eq/P))
    c=[sp.simplify(eq.coeff(sp.Derivative(h,(u,k))) if k>0 else eq.subs({sp.Derivative(h,(u,2)):0,sp.Derivative(h,u):0}).coeff(h)) for k in range(3)]
    # common factor to make finite at u=1
    c=[sp.factor(ci) for ci in c]
    return c
import pickle
c=build(-2,sp.Rational(1,10**6),0,2,4)
for k,ci in enumerate(c): print(k, sp.limit(sp.simplify(ci*(1-u)**4),u,1) if True else '')
for p in (0,1,2,3):
    print(p,[sp.limit(sp.simplify(ci*(1-u)**p),u,1) for ci in c])

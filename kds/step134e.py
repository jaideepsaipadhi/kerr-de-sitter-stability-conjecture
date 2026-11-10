import sympy as sp, numpy as np
z,t,z0,A,B,v,be=sp.symbols('z t z0 A B v beta')
s=z0-z; mu=z*(z-1)*(t-z); L=A/z+B-3*z; F=mu/s**2
nu,dn=sp.fraction(sp.factor(sp.together(sp.diff(mu*sp.diff(F,z),z)-2*s*sp.diff(L*F,z))))
x0=z0-1; N=sp.expand(nu.subs(z,z0-x0*v/(1+v))*(1+v)**6)
num=sp.Poly(sp.expand(sp.cancel(N)),v)
T=[sp.expand(num.coeff_monomial(v**k).subs(B,3*z0-A/z0+be)*z0) for k in range(7)]
Ac2=z0*((t-1)/(2*x0**2)-3)
f=[sp.lambdify((t,z0,A,be),Tk) for Tk in T]
rng=np.random.default_rng(3); worst={k:-1e99 for k in (4,5,6)}; bad=0
for _ in range(200000):
    tt=1+10**rng.uniform(-4,5); zz=1+(tt-1)*rng.uniform(1e-4,1-1e-4)
    a0=max(0.0,zz*((tt-1)/(2*(zz-1)**2)-3))
    for k in (4,5,6):
        val=f[k](tt,zz,a0,0.0); 
        if val>0: bad+=1; 
        worst[k]=max(worst[k],val/abs(f[k](tt,zz,a0+1,1)-val+1e-300))
print('bad',bad,worst)
import pickle; pickle.dump(T,open('T134.pkl','wb'))

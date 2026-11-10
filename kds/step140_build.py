import sympy as sp, pickle, time
t0=time.time()
z,t,z0,A,B,S,be,sg=sp.symbols('z t z0 A B S beta sigma',positive=True)
s=z0-z; mu=z*(z-1)*(t-z); L=A/z+B-3*z; F=mu/s**2
nu,dn=sp.fraction(sp.factor(sp.together(sp.diff(mu*sp.diff(F,z),z)-2*s*sp.diff(L*F,z))))
N=sp.expand(nu.subs(z,z0-S).subs(B,3*z0-A/z0+be)*z0)
cb=N.coeff(be); N0=sp.expand(N-cb*be)
beta_star=sp.factor(-N0/cb)
print('beta* =',sp.count_ops(beta_star),'ops; cb=',sp.factor(cb), round(time.time()-t0,1))
# barrier G/s at z = z0 - S*sigma, with w0=H/S at zhat=z0-S
zh=z0-S; muh=zh*(zh-1)*(t-zh); H=2*muh+S*sp.diff(mu,z).subs(z,zh); w0=H/S
q=z0-S*sg; sq=S*sg; muq=q*(q-1)*(t-q); Bv=3*z0-A/z0+be
Lq=A/q+Bv-3*q; Lp=A/q**2+3
G=-2*muq*w0/S**2+w0**2*sq**3/S**4-2*Lq*w0*sq**2/S**2+2*muq*Lp
Gn,Gd=sp.fraction(sp.together(G))
print('G den',sp.factor(Gd), round(time.time()-t0,1))
pickle.dump((beta_star,cb,N0,Gn,Gd),open('P140.pkl','wb'))
print('G num deg in beta',sp.Poly(Gn,be).degree(), 'total terms',len(sp.expand(Gn).args), round(time.time()-t0,1))

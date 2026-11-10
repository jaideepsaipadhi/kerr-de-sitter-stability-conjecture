import sympy as sp, pickle
z,t,z0,A,B=sp.symbols('z t z0 A B')
s=z0-z; mu=z*(z-1)*(t-z); L=A/z+B-3*z; F=mu/s**2
nu,dn=sp.fraction(sp.factor(sp.together(sp.diff(mu*sp.diff(F,z),z)-2*s*sp.diff(L*F,z))))
N=sp.expand(nu); H=2*mu+s*sp.diff(mu,z); T=sp.expand(H-s**2*L)
NB=N.coeff(B); TB=T.coeff(B)
# Phi = T*NB - TB*N  (B cancels); with NB<0: T - kappa N >=0 where kappa=TB/NB>0  <=>  T*NB - TB*N <= 0
Phi=sp.factor(sp.expand(T*NB-TB*N))
print('B in Phi?',Phi.has(B)); print(Phi)
pickle.dump(Phi,open('P156.pkl','wb'))

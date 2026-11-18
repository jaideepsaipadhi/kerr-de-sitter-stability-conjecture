# C3 builder: coefficients t_0..t_6 of the sextic (Lemma "One sign change") and the
# vertex values g1, g2 of the triangle argument. Writes T134.pkl and G135.pkl,
# which step135g.py then certifies.
import sympy as sp, pickle
z,t,z0,A,B,v,be=sp.symbols('z t z0 A B v beta')
s=z0-z; mu=z*(z-1)*(t-z); L=A/z+B-3*z; F=mu/s**2
nu,dn=sp.fraction(sp.factor(sp.together(sp.diff(mu*sp.diff(F,z),z)-2*s*sp.diff(L*F,z))))
x0=z0-1; N=sp.expand(nu.subs(z,z0-x0*v/(1+v))*(1+v)**6)
num=sp.Poly(sp.expand(sp.cancel(N)),v)
T=[sp.expand(num.coeff_monomial(v**k).subs(B,3*z0-A/z0+be)*z0) for k in range(7)]
pickle.dump(T,open('T134.pkl','wb'))
a0p=z0*((t-1)/(2*x0**2)-3)
G={}
for case,a0 in (('I',a0p),('II',sp.Integer(0))):
  for j in (2,3,4):
    Tj,Tk=T[j],T[j+1]
    a=sp.together(Tj.subs({A:a0,be:0})); b=sp.together(Tk.subs({A:a0,be:0}))
    pA,pB=Tj.coeff(A),Tj.coeff(be); qA,qB=Tk.coeff(A),Tk.coeff(be)
    # g1 = |qA| * t_j at (a0 + b/|qA|, 0), g2 = |qB| * t_j at (a0, b/|qB|)
    G[(case,j)]=(a,b,sp.together(-a*qA+pA*b),sp.together(-a*qB+pB*b),pA,pB,qA,qB)
pickle.dump(G,open('G135.pkl','wb'))
print('built')

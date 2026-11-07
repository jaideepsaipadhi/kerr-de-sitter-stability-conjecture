# Principal symbol (fixed mode, variables (r,theta)) of the DGS coupled system, eqs (4),(7).
import sympy as sp
X,Y,p,q,D=sp.symbols('X Y p q Delta',positive=True)  # X=sqrt(Delta)xi_r, Y=xi_theta, p=Q^2 a_- rbar*, q=Q^2 a_+ rbar*
Mt=sp.Matrix([[(1+p)*X**2+(1-q)*Y**2, -(p+q)*sp.sqrt(D)*X*Y],[-(p+q)*X*Y/sp.sqrt(D),(1-q)*X**2+(1+p)*Y**2]])
print("det =",sp.factor(Mt.det()))
print("eigenvalues =",[sp.factor(e) for e in Mt.eigenvals()])
r,a,th,M,Q=sp.symbols('r a theta M Q',positive=True); rb=r+sp.I*a*sp.cos(th); rbs=r-sp.I*a*sp.cos(th)
am=1/(3*(rb**2*M-rb*Q**2)-Q**2*rbs); ap=1/(3*(rb**2*M-rb*Q**2)+Q**2*rbs)
print("1+p =",sp.simplify(1+Q**2*am*rbs)," ; 1-q =",sp.simplify(1-Q**2*ap*rbs))

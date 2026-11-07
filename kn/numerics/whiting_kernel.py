# Step 1: Euler-type kernel for the confluent Heun form of the radial Teukolsky equation.
# L_z = z(z-1)d^2 + (g(z-1)+d z + e z(z-1)) d + (al z - q). Seek K=exp(p x z) with
# L_z^dagger K = M_x K, M_x confluent Heun in x with q' = q + c (c independent of q => lambda-free).
import sympy as sp
x,z,p,g,d,e,al,q=sp.symbols('x z p gamma delta epsilon alpha q')
K=sp.exp(p*x*z)
A=z*(z-1); B=g*(z-1)+d*z+e*z*(z-1); C=al*z-q
Ldag=sp.diff(A*K,z,2)-sp.diff(B*K,z)+C*K
g2,d2,e2,al2,c=sp.symbols('g2 d2 e2 al2 c')
M=x*(x-1)*sp.diff(K,x,2)+(g2*(x-1)+d2*x+e2*x*(x-1))*sp.diff(K,x)+(al2*x-(q+c))*K
expr=sp.expand(sp.simplify((Ldag-M)/K))
eqs=sp.Poly(expr,x,z).coeffs()
for s in sp.solve(eqs,[p,g2,d2,e2,al2,c],dict=True): print(s)

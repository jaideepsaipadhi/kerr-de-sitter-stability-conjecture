# C6 (exact): divergence of the self-dual 2-form phi*(e0^e1 - i e2^e3), phi=(r+ip)/Omega, on the
# accelerating Plebanski-Demianski metric with NUT parameter 0. Prints xi = (3, -3 i alpha, 0, 0) in (tau,sigma,r,p),
# exactly, for all m, k, epsilon, Lambda, alpha.
import sympy as sp
tau,sig,r,p=sp.symbols('tau sigma r p',real=True); X=[tau,sig,r,p]
al,m,k,e,L=sp.symbols('alpha m k epsilon Lambda',real=True); n=0
P=k+2*n*p-e*p**2+2*al*m*p**3-(al**2*k+L/3)*p**4
Q=k-2*m*r+e*r**2-2*al*n*r**3-(al**2*k+L/3)*r**4
S=r**2+p**2; Om=1-al*p*r
w1=[1,-p**2,0,0]; w2=[1,r**2,0,0]
g=sp.zeros(4)
for i in range(4):
    for j in range(4):
        g[i,j]=(Q/S*w1[i]*w1[j]-P/S*w2[i]*w2[j])/Om**2
g[2,2]=-S/Q/Om**2; g[3,3]=-S/P/Om**2
gi=sp.simplify(g.inv()); sq=sp.simplify(sp.sqrt(sp.factor(-g.det())))
print('sqrtg',sq)
def wedge(a,b): return sp.Matrix(4,4,lambda i,j:a[i]*b[j]-a[j]*b[i])
e01=wedge([1,-p**2,0,0],[0,0,1,0])/Om**2       # (dtau - p^2 dsig) ^ dr
e23=wedge([0,0,0,1],[1,r**2,0,0])/Om**2        # dp ^ (dtau + r^2 dsig)
for sgn in (1,):
  for phs in (-1,):
    phi=(r+sgn*sp.I*p)/Om
    F=phi*(e01+phs*sp.I*e23)
    Fu=gi*F*gi.T
    xi=[sp.simplify(sum(sp.diff(sq*Fu[b,a],X[b]) for b in range(4))/sq) for a in range(4)]
    print(sgn,phs,[sp.factor(x) for x in xi])

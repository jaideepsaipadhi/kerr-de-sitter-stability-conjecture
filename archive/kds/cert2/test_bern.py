# Soundness test of bern_core: Bernstein bounds after random subdivisions must enclose the true range.
import sympy as sp, random, numpy as np, itertools
from bern_core import to_bern, split, mn, mx
x,y,z=sp.symbols('x y z'); V=(x,y,z); random.seed(7); bad=0
for trial in range(40):
    P=sp.Poly(sum(random.randint(-9,9)*x**i*y**j*z**k for i in range(4) for j in range(3) for k in range(3) if random.random()<0.5),*V)
    if P.is_zero: continue
    B=to_bern(P,V); box=[(0.0,1.0)]*3
    for _ in range(6):
        k=random.randrange(3); side=random.randrange(2); B=split(B,k)[side]
        a,b=box[k]; m=(a+b)/2; box[k]=(a,m) if side==0 else (m,b)
    # Bernstein arrays are scaled by positive integers; compare SIGN-relevant bounds via normalisation by corner value
    f=sp.lambdify(V,P.as_expr(),'math')
    pts=[tuple(random.uniform(*box[d]) for d in range(3)) for _ in range(400)]+list(itertools.product(*box))
    vals=[f(*p) for p in pts]
    # scale factor: compare corner value with corner Bernstein coefficient (exact relation: coefficient = scale * value)
    c0=B[(0,0,0)]; v0=f(box[0][0],box[1][0],box[2][0])
    if v0==0: continue
    sc=float(c0)/v0
    lo,hi=float(mn(B))/sc,float(mx(B))/sc
    if sc<0: lo,hi=hi,lo
    if min(vals)<lo-1e-9*abs(lo)-1e-12 or max(vals)>hi+1e-9*abs(hi)+1e-12 or sc<=0: bad+=1; print('FAIL',trial,lo,hi,min(vals),max(vals),sc)
print('enclosure failures:',bad)

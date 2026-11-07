import sympy as sp, mpmath as mp
exec(open('nh_coupled.py').read().split('def lead')[0].replace("f2,f1=sp.Function('f2')(th),sp.Function('f1')(th)",
 "f2=sp.exp(sp.cos(th))*sp.sin(th)**3; f1=sp.exp(2*sp.cos(th))*sp.sin(th)**2"))
mp.mp.dps=40; hv=sp.Rational(3,2)+sp.I*sp.Rational(7,10)
for d in [-2,-1,0,1,2]:
    p2=x**h*f2; p1=x**(h+d)*f1
    res=[]
    for k,e in enumerate(ops(p2,p1)):
        base=h if k<3 else h+d
        g=sp.lambdify((x,),(e/x**base).subs({a:sp.Rational(7,10),m:2,h:hv,th:sp.Rational(9,10)}),'mpmath')
        v1,v2=abs(g(mp.mpf('1e-6'))),abs(g(mp.mpf('1e-8')))
        res.append(round(float(mp.log(v2/v1)/mp.log(mp.mpf('1e-2'))),2) if v1>0 else None)
    print(f"d={d}: psi-2 eq powers(F,G,H)={res[:3]}   psi-1 eq powers(F,G,H)={res[3:]}")

import sympy as sp, random
r,M,a,Q,W,K=sp.symbols('r M a Q Omega K',real=True)
D=r**2-2*M*r+a**2+Q**2; P=(r**2+a**2)*W-a
G=2*sp.diff(P,r)*D-sp.diff(D,r)*P
h=r**3-3*M*r**2+(a**2+2*Q**2)*r+M*a**2
print("G/2 identity:",sp.simplify(G/2-(W*h+a*(r-M)))==0)
rp=sp.symbols('r_p',positive=True)
sub={Q**2:2*M*rp-rp**2-a**2}
k=r**3-3*M*r**2+(a**2+Q**2+2*M*rp)*r+M*(a**2+Q**2)-2*M**2*rp
print("k identity:",sp.simplify((h+(2*M*rp-Q**2)*(r-M))-k)==0)
print("k(r+)=0:",sp.simplify(k.subs(r,rp).subs(sub))==0)
print("k'(M):",sp.factor(sp.diff(k,r).subs(r,M).subs(sub)))
R=P**2-D*K; print("R'=0 reduction:",sp.simplify(sp.diff(R,r).subs(K,P**2/D)*D-P*G)==0)
bad=0
for _ in range(20000):
    aa=random.uniform(1e-3,.999);qq=random.uniform(0,(1-aa*aa)**.5*.999);rpp=1+(1-aa*aa-qq*qq)**.5
    w=random.uniform(0,1)*aa/(rpp**2+aa**2);x=rpp+random.expovariate(0.3)+1e-9
    if float(G.subs({M:1,a:aa,Q:qq,W:w,r:x}))<=0: bad+=1
print("sampled violations:",bad)

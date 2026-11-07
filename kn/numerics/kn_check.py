import sympy as sp, random
r,M,a,W,Q=sp.symbols('r M a Omega Q',positive=True)
D=r**2-2*M*r+a**2+Q**2
P=(r**2+a**2)*W-a
G=sp.expand(2*sp.diff(P,r)*D-P*sp.diff(D,r))   # F'=P*G/D^2, F=P^2/D
rp=M+sp.sqrt(M**2-a**2-Q**2); wH=a/(rp**2+a**2)
# equatorial: is {F = (1-a W)^2} the same as {V=dt+W dphi null at theta=pi/2}?
S=r**2  # Sigma at equator
gtt=-(1-2*M*r/S+Q**2/S); gtp=-a*(2*M*r-Q**2)/S
gpp=(r**2+a**2+a**2*(2*M*r-Q**2)/S)
VV=sp.simplify((gtt+2*W*gtp+W**2*gpp)*r**2)
print("V.V*r^2 + (P^2-(1-aW)^2 D) =",sp.simplify(VV+(P**2-(1-a*W)**2*D)))
# numeric scan: G>0 on (r+,inf) for W in (0,wH]?
bad=0
for _ in range(4000):
    m=1.0; aa=random.uniform(0,.999); qq=random.uniform(0,(1-aa**2)**.5*.999)
    rpp=1+(1-aa**2-qq**2)**.5; w=random.uniform(0,1)*aa/(rpp**2+aa**2)
    g=sp.lambdify(r,G.subs({M:1,a:aa,Q:qq,W:w}))
    for x in [rpp*(1+1e-6)+k*0.05 for k in range(400)]:
        if g(x)<=0: bad+=1;break
print("violations of G>0 (Kerr-Newman, Omega in superradiant band):",bad)
print("V.V*r^2 - (P^2-(1-aW)^2 D) =",sp.simplify(VV-(P**2-(1-a*W)**2*D)))
print("ratio:",sp.factor(sp.simplify(VV/(P**2-(1-a*W)**2*D))))

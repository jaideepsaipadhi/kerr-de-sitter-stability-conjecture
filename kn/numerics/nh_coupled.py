# Leading near-horizon balance of the DGS coupled system (their eqs (4),(7),(8)) on EXTREMAL KN, omega = m*Omega_H.
import sympy as sp
x,th,h,a,m=sp.symbols('x theta h a m'); Q=sp.sqrt(1-a**2); M=1
r=1+x; Dl=x**2; w=m*a/(1+a**2)
Kr=a*m-(r**2+a**2)*w; Kt=m/sp.sin(th)-a*w*sp.sin(th)
rb=r+sp.I*a*sp.cos(th); rbs=r-sp.I*a*sp.cos(th)
al=lambda sg: 1/(3*(rb**2*M-rb*Q**2)+sg*Q**2*rbs)
D =lambda j,u: sp.diff(u,x)+sp.I*Kr/Dl*u+2*j*(r-M)/Dl*u
Dd=lambda j,u: sp.diff(u,x)-sp.I*Kr/Dl*u+2*j*(r-M)/Dl*u
L =lambda j,u: sp.diff(u,th)+Kt*u+j*sp.cot(th)*u
Ld=lambda j,u: sp.diff(u,th)-Kt*u+j*sp.cot(th)*u
iasin=sp.I*a*sp.sin(th)
f2,f1=sp.Function('f2')(th),sp.Function('f1')(th)
def ops(p2,p1):
    p2=p2; p1=p1
    F2=Dl*Dd(-1,D(0,p2))+L(-1,Ld(2,p2))-6*sp.I*w*rb*p2
    G2=Dl*Dd(-1,al(-1)*rbs*D(0,p2))-3*Dl*Dd(-1,al(-1)*p2)-L(-1,al(1)*rbs*Ld(2,p2))+3*L(-1,al(1)*iasin*p2)
    H2=-Dl*Dd(-1,al(-1)*rbs*L(-1,p1))-3*Dl*Dd(-1,al(-1)*iasin*p1)-L(-1,al(1)*rbs*Dl*Dd(-1,p1))-3*L(-1,al(1)*Dl*p1)
    F1=Dl*D(1,Dd(-1,p1))+Ld(2,L(-1,p1))-6*sp.I*w*rb*p1
    G1=-D(0,al(1)*rbs*Dl*Dd(-1,p1))-3*D(0,al(1)*Dl*p1)+Ld(2,al(-1)*rbs*L(-1,p1))+3*Ld(2,al(-1)*iasin*p1)
    H1=-D(0,al(1)*rbs*Ld(2,p2))+3*D(0,al(1)*iasin*p2)-Ld(2,al(-1)*rbs*D(0,p2))+3*Ld(2,al(-1)*p2)
    return F2,G2,H2,F1,G1,H1
def lead(e,p0):
    e=sp.simplify(e/x**p0)
    s=sp.series(e,x,0,1).removeO()
    for k in range(-4,1):
        c=sp.simplify(s.coeff(x,k))
        if c!=0: return k
    return None
av=sp.Rational(7,10); mv=2
for d in [-2,-1,0,1,2]:
    p2=x**h*f2; p1=x**(h+d)*f1
    F2,G2,H2,F1,G1,H1=[sp.expand(o.subs({a:av,m:mv})) for o in ops(p2,p1)]
    o=[lead(e,h) for e in (F2,G2,H2)]+[lead(e,h+d) for e in (F1,G1,H1)]
    print(f"d=h1-h2={d}: psi-2 eq leading x-powers F,G,H = {o[:3]};  psi-1 eq F,G,H (rel. x^h1) = {o[3:]}")

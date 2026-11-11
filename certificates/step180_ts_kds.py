# Exact check: on KdS, psi = Delta^2 (D^+)^4 (Delta^2 R) solves the s=-2 radial equation whenever R solves the s=+2
# equation with the same lambda. Radial eq (step1_radial.py): Delta R'' + (s+1) Delta' R' + Q_s R = 0,
# Q_s = (Xi^2 K^2 - i s Xi K Delta')/Delta + 4 i s Xi sg r - (2/3) Lam (s+1)(2s+1) r^2 + s(1-alpha) - lam,
# K = sg (r^2+a^2) - a m,  D^+ = d/dr + i Xi K/Delta.  Fully symbolic in (Lam, a, M, sg, m, lam).
# Elements: (A, B, k) meaning (A R + B R')/Delta^k with A, B polynomials over QQ(i).
import sympy as sp, time
t0=time.time()
r,L,a,M,sg,m,lam=sp.symbols('r Lam a M sg m lam')
G=(r,L,a,M,sg,m,lam); P=lambda e: sp.Poly(e,*G,domain='QQ_I')
I=sp.I; al=L*a**2/3; Xi=1+al
Dl=(r**2+a**2)*(1-L*r**2/3)-2*M*r; Dlp=sp.diff(Dl,r); K=sg*(r**2+a**2)-a*m
def Qn(s):  # Delta*Q_s, polynomial
    return sp.expand(Xi**2*K**2-I*s*Xi*K*Dlp+Dl*(4*I*s*Xi*sg*r-sp.Rational(2,3)*L*(s+1)*(2*s+1)*r**2+s*(1-al)-lam))
D_,Dp_,Q2=P(Dl),P(Dlp),P(Qn(2)); XK=P(I*Xi*K)
def dr(p): return p.diff(r)
def deriv(e):   # d/dr, reducing R'' via the s=+2 equation: R'' = -(3 Delta' R' + Q R)/Delta, Q = Q2/Delta
    A,B,k=e
    A2=dr(A)*D_**2 - k*Dp_*D_*A - B*Q2
    B2=A*D_**2 + dr(B)*D_**2 - 3*Dp_*D_*B - k*Dp_*D_*B
    return (A2,B2,k+2)
def Dplus(e):
    A,B,k=deriv(e); A0,B0,k0=e   # add i Xi K/Delta * e  -> over Delta^{k0+1}; bring to Delta^k (k = k0+2)
    return (A+XK*A0*D_, B+XK*B0*D_, k)
e=(D_**2, P(0), 0)
for j in range(4): e=Dplus(e); print('D+',j,round(time.time()-t0,1),flush=True)
e=(D_**2*e[0], D_**2*e[1], e[2])
e1=deriv(e); e2=deriv(e1)
# s=-2 operator: Delta psi'' - Delta' psi' + Q_{-2} psi, with Q_{-2} = Qm/Delta. Common denominator Delta^{k+4}... 
k0,k1,k2=e[2],e1[2],e2[2]; Qm=P(Qn(-2)); N=k2+1
def lift(x,kk,mult): return x*mult*D_**(N-kk)   # x*mult/Delta^kk over Delta^N
RA=lift(e2[0],k2,D_)-lift(e1[0],k1,Dp_)+lift(e[0],k0+1,Qm)
RB=lift(e2[1],k2,D_)-lift(e1[1],k1,Dp_)+lift(e[1],k0+1,Qm)
print('R coefficient zero:',RA.is_zero,' R\' coefficient zero:',RB.is_zero, round(time.time()-t0,1))

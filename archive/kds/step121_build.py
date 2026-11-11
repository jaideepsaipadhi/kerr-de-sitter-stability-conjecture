# (P2) certificate polynomials in (T = 1/t, ze, sg = S/x0):
#   NA0 := N(S; A = 0, B = floor)  [case A_K <= 0],  NAK := N(S; A = A_K(S), B = floor) [case A_K >= 0],
#   AKn := numerator of A_K (case gate), BC := B-coefficient of N (must be < 0).
import sympy as sp, pickle, time
T,ze,sg=sp.symbols('T ze sg',positive=True); V=(T,ze,sg); t0=time.time()
t=1/T; z0=1+ze*(t-1); x0=z0-1; S=sg*x0; z=z0-S
mu=z*(z-1)*(t-z); mup=-3*z**2+2*(1+t)*z-t; H=2*mu+S*mup
mu0=z0*x0*(t-z0); mup0=-3*z0**2+2*(1+t)*z0-t; Bf=sp.Rational(11,2)*t-sp.Rational(5,4)
def N(A,B):
    c=[6*mu0**2,-6*mu0*mup0,4*(z0-1)*(z0-t)*(A+B*z0)+6*t**2*z0**2-6*t**2*z0+t**2-8*t*z0**3+12*t*z0**2-6*t*z0+3*z0**4-8*z0**3+6*z0**2,
       2*(A*(t+1-2*z0)+B*mup0+3*z0*(-3*t*z0+2*t+4*z0**2-3*z0)),0,2*(B+4*t-15*z0+4),15]
    return sum(c[k]*S**k for k in range(7))
AK=(z**2*H-3*z**2*S**3)/S**3
out={}
def rec(nm,ex):
    nu,dn=sp.fraction(sp.factor(sp.together(ex))); out[nm]=(sp.Poly(sp.expand(nu),*V),dn)
    print(nm,[out[nm][0].degree(v) for v in V],'den',str(dn)[:70],round(time.time()-t0,1),flush=True)
rec('AK',AK); rec('NA0',N(0,Bf)); rec('NAK',N(AK,Bf)); rec('BC',4*z0*x0*(z0-t)+2*mup0*S+2*S**3)
pickle.dump(out,open('/home/claude/kds/p2_3d.pkl','wb'))

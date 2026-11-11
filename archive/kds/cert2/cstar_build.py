import sympy as sp, pickle, time
from ratpoly import R
T,ze,y=sp.symbols('T ze y',positive=True); V=(T,ze,y); t0=time.time()
P=lambda e: R(sp.Poly(e,*V))
Tp=P(T); one=P(1)
t=one/Tp; x0=P(ze)*(t-1); z0=x0+1; x=P(y)*x0; z=x+1; s=(1-P(y))*x0; mu0=z0*x0*(t-z0)
mu=z*x*(t-z); mup=z*z*(-3)+(t+1)*z*2-t
v=1-P(y); q=v*sp.Rational(9,20)/(x*z); qp=((one/x0)/(x*z)*(-1)-v*(z*2-1)/((x*z)*(x*z)))*sp.Rational(9,20)
L0s=mu0*sp.Rational(5,4)/(s*s); kap=s/z; L1=(z0/z)*(L0s+z0*3)-z*3
g0=(one+one/kap)/z; C=(L1/kap-z*6)/z; fl=t*sp.Rational(11,2)-sp.Rational(5,4)
E=s*(g0*sp.Rational(5455,10000)+q); Lhi=L1-kap*fl
Bhi=L0s+z0*3; Bpp=(z0*(L0s+z0*3)-3-(t-1)/(x0*2))/(z0-1); Bc2=L0s+z0*3-((t-1)/(x0*x0*2)-3)
out={}
def rec(ex,nm):
    out[nm]=(ex.n,ex.d); print(nm,'num deg',[ex.n.degree(v_) for v_ in V],'den deg',[ex.d.degree(v_) for v_ in V],round(time.time()-t0,1),flush=True)
for br,c in (('pos',sp.Rational(17,10)),('neg',sp.Rational(17,11))):
    G1c=mup*c*sp.Rational(1,2)+mu*q*sp.Rational(17,10)
    a2=mu*C*C*sp.Rational(2295,1000); a1=-(C*(mu*g0*sp.Rational(459,100)+G1c))
    rec(a2*a2*E*27+a1**3,'P_'+br)
rec(E*Lhi**3-mu*C*C*sp.Rational(2295,1000),'region')
rec(Bhi-fl,'f1'); rec(Bpp-fl,'f2'); rec(Bc2-fl,'f3'); rec(Lhi,'Lpos'); rec(mup,'mup')
pickle.dump(out,open('cstar.pkl','wb'))

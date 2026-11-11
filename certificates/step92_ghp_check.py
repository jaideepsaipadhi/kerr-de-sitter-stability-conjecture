# Item 92: machine check (KdS, Carter/Kinnersley tetrad, numeric derivatives at 40 digits) of the GHP relations of item 91.
# Conventions fixed by this check: Th' = Delta - p gamma - q gammabar, eth' = deltabar - p alpha - q betabar.
# Item 92: machine check of the type D GHP relations on KdS (Carter tetrad), numeric point evaluation.
import sympy as sp, time
t0=time.time()
t,r,th,ph=sp.symbols('t r th ph',real=True); X=[t,r,th,ph]
M,a,L=sp.Rational(1),sp.Rational(3,5),sp.Rational(1,20)
al=L*a**2/3; Xi=1+al
Dr=(r**2+a**2)*(1-L*r**2/3)-2*M*r; Dth=1+al*sp.cos(th)**2
rho2=r**2+a**2*sp.cos(th)**2
s=sp.sin(th)
# signature (+,-,-,-): g = Dr/(rho2 Xi^2)(dt - a s^2 dph)^2 - Dth s^2/(rho2 Xi^2)(a dt-(r^2+a^2)dph)^2 - rho2(dr^2/Dr + dth^2/Dth)
e1=[1,0,0,-a*s**2]; e2=[a,0,0,-(r**2+a**2)]
g=sp.zeros(4)
for i in range(4):
    for j in range(4):
        g[i,j]=Dr/(rho2*Xi**2)*e1[i]*e1[j]-Dth*s**2/(rho2*Xi**2)*e2[i]*e2[j]
g[1,1]=-rho2/Dr; g[2,2]=-rho2/Dth
gi=g.inv()
pt={r:sp.Rational(3),th:sp.Rational(7,10)}
c=r-sp.I*a*sp.cos(th)
l=[Xi*(r**2+a**2)/Dr,1,0,Xi*a/Dr]                 # Kinnersley-type (l.n = 1)
n=[Xi*(r**2+a**2)/(2*rho2),-Dr/(2*rho2),0,Xi*a/(2*rho2)]
m=[sp.I*Xi*a*s/(sp.sqrt(2*Dth)*(r+sp.I*a*sp.cos(th))),0,sp.sqrt(Dth)/(sp.sqrt(2)*(r+sp.I*a*sp.cos(th))),sp.I*Xi/(sp.sqrt(2*Dth)*s*(r+sp.I*a*sp.cos(th)))]
mb=[sp.conjugate(x) for x in m]
dot=lambda u,v: sum(g[i,j]*u[i]*v[j] for i in range(4) for j in range(4))
print('norms', [sp.N(sp.simplify(dot(u,v)).subs(pt)) for u,v in [(l,n),(m,mb),(l,l),(m,m),(l,m)]], flush=True)
Gam=[[[sum(gi[k,q]*(sp.diff(g[q,i],X[j])+sp.diff(g[q,j],X[i])-sp.diff(g[i,j],X[q]))/2 for q in range(4)) for j in range(4)] for i in range(4)] for k in range(4)]
low=lambda u:[sum(g[i,j]*u[j] for j in range(4)) for i in range(4)]
def cov(u):  # nabla_b u_a as matrix [a][b]
    ul=low(u)
    return [[sp.diff(ul[i],X[j])-sum(Gam[k][i][j]*ul[k] for k in range(4)) for j in range(4)] for i in range(4)]
Nl,Nn,Nm,Nmb=cov(l),cov(n),cov(m),cov(mb)
c3=lambda A,u,v: sum(A[i][j]*u[i]*v[j] for i in range(4) for j in range(4))  # u^a v^b nabla_b A_a
kap=c3(Nl,m,l); rho=c3(Nl,m,mb); sig=c3(Nl,m,m); tau=c3(Nl,m,n)
nu=-c3(Nn,mb,n); pi=-c3(Nn,mb,l); lam=-c3(Nn,mb,mb); mu=-c3(Nn,mb,m)
eps=(c3(Nl,n,l)+c3(Nmb,m,l))/2; beta=(c3(Nl,n,m)+c3(Nmb,m,m))/2
alp=(c3(Nl,n,mb)+c3(Nmb,m,mb))/2; gam=(c3(Nl,n,n)+c3(Nmb,m,n))/2

import mpmath as mp
mp.mp.dps=40
names=['kap','rho','sig','tau','nu','pi','lam','mu','eps','beta','alp','gam']
exprs=[kap,rho,sig,tau,nu,pi,lam,mu,eps,beta,alp,gam]
fn={k:sp.lambdify((r,th),e,'mpmath') for k,e in zip(names,exprs)}
vec={k:[sp.lambdify((r,th),x,'mpmath') for x in u] for k,u in [('l',l),('n',n),('m',m)]}
print('lambdified',round(time.time()-t0,1),flush=True)
R0,T0=mp.mpf(3),mp.mpf('0.7')
SC={k:f(R0,T0) for k,f in fn.items()}
def V(F,u):  # directional derivative of scalar function F(r,th) along vector u at point (only r, th components act)
    ur=vec[u][1](R0,T0); ut=vec[u][2](R0,T0)
    return ur*mp.diff(lambda x:F(x,T0),R0)+ut*mp.diff(lambda y:F(R0,y),T0)
def conjv(u): return u
def Vb(F):  # mbar direction
    ur=mp.conj(vec['m'][1](R0,T0)); ut=mp.conj(vec['m'][2](R0,T0))
    return ur*mp.diff(lambda x:F(x,T0),R0)+ut*mp.diff(lambda y:F(R0,y),T0)
c=mp.conj
e,b,al_,ga=SC['eps'],SC['beta'],SC['alp'],SC['gam']
def Th(F,p,q):  return V(F,'l')-(p*e+q*c(e))*F(R0,T0)
SG=[1,1]
def Thp(F,p,q): return V(F,'n')+-(p*ga+q*c(ga))*F(R0,T0)
def Ed(F,p,q):  return V(F,'m')-(p*b+q*c(al_))*F(R0,T0)
def Edp(F,p,q): return Vb(F)+(-p*al_-q*c(b))*F(R0,T0)
RHO=fn['rho']; TAU=fn['tau']; RP=lambda x,y:-fn['mu'](x,y); TP=lambda x,y:-fn['pi'](x,y)
rho_,tau_,rp_,tp_=RHO(R0,T0),TAU(R0,T0),RP(R0,T0),TP(R0,T0)
if True:
 if True:
  chk={
   'Th rho - rho^2': Th(RHO,1,1)-rho_**2,
   'Ed rho - (rho-rhob)tau': Ed(RHO,1,1)-(rho_-c(rho_))*tau_,
   'Th tau - (tau - taubp)rho': Th(TAU,1,-1)-(tau_-c(tp_))*rho_,
   'Ed tau - tau^2': Ed(TAU,1,-1)-tau_**2,
   'Thp rp - rp^2': Thp(RP,-1,-1)-rp_**2,
   '[1] Th rp - Thp rho - (tau taub - tp tpb)': Th(RP,-1,-1)-Thp(RHO,1,1)-(tau_*c(tau_)-tp_*c(tp_)),
   '[2] Ed tp - Edp tau - (rho rpb - rhob rp)': Ed(TP,-1,1)-Edp(TAU,1,-1)-(rho_*c(rp_)-c(rho_)*rp_),
   '[3] Th tp - Edp rho': Th(TP,-1,1)-Edp(RHO,1,1),
   '[4] Thp tau - Ed rp': Thp(TAU,1,-1)-Ed(RP,-1,-1),
  }
  for k,v in chk.items(): print('   ',k, mp.nstr(abs(v),5), flush=True)

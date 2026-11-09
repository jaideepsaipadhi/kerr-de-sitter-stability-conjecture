# Hostile check: the GHP tetrad-derivative table (item 93(b)) on the accelerating NUT Lambda-PD metric,
# where tau, rho etc. are generic complex numbers (on KdS tau is imaginary, so tau and -taubar coincide).
import sympy as sp, mpmath as mp, time
mp.mp.dps=40
tau_,sig_,r,p=sp.symbols('tau sig r p',real=True); X=[tau_,sig_,r,p]
Lam=sp.Rational(1,20); a=sp.Rational(3,5); m_,n_,al_=1,sp.Rational(1,10),sp.Rational(1,20)
e_=1-Lam*a**2/3; k_=a**2
P=k_+2*n_*p-e_*p**2+2*al_*m_*p**3-(al_**2*k_+Lam/3)*p**4
Q=k_-2*m_*r+e_*r**2-2*al_*n_*r**3-(al_**2*k_+Lam/3)*r**4
S=r**2+p**2; Om=1-al_*p*r
w1=[1,-p**2,0,0]; w2=[1,r**2,0,0]
g=sp.zeros(4)
for i in range(4):
    for j in range(4): g[i,j]=(Q/S*w1[i]*w1[j]-P/S*w2[i]*w2[j])/Om**2
g[2,2]=-S/Q/Om**2; g[3,3]=-S/P/Om**2
gi=g.inv(); V1=[r**2/S,-1/S,0,0]; V2=[p**2/S,1/S,0,0]; A=Om*sp.sqrt(S/(2*Q))
l=[A*V1[i]+A*(Q/S)*(1 if i==2 else 0) for i in range(4)]; n=[A*V1[i]-A*(Q/S)*(1 if i==2 else 0) for i in range(4)]
m=[Om/sp.sqrt(2)*(sp.sqrt(P/S)*(1 if i==3 else 0)+sp.I*sp.sqrt(S/P)*V2[i]) for i in range(4)]; mb=[sp.conjugate(x) for x in m]
Gam=[[[sum(gi[k,q]*(sp.diff(g[q,i],X[j])+sp.diff(g[q,j],X[i])-sp.diff(g[i,j],X[q])) for q in range(4))/2 for j in range(4)] for i in range(4)] for k in range(4)]
low=lambda u:[sum(g[i,j]*u[j] for j in range(4)) for i in range(4)]
def cov(u):
    ul=low(u); return [[sp.diff(ul[i],X[j])-sum(Gam[k][i][j]*ul[k] for k in range(4)) for j in range(4)] for i in range(4)]
Nl,Nn,Nmb=cov(l),cov(n),cov(mb)
c3=lambda Aa,u,v: sum(Aa[i][j]*u[i]*v[j] for i in range(4) for j in range(4))
E={'rho':c3(Nl,m,mb),'tau':c3(Nl,m,n),'mu':-c3(Nn,mb,m),'pi':-c3(Nn,mb,l),'eps':(c3(Nl,n,l)+c3(Nmb,m,l))/2,
   'beta':(c3(Nl,n,m)+c3(Nmb,m,m))/2,'alp':(c3(Nl,n,mb)+c3(Nmb,m,mb))/2,'gam':(c3(Nl,n,n)+c3(Nmb,m,n))/2}
R0,P0=mp.mpf(3),mp.mpf('0.36')
sc={k:sp.lambdify((r,p),v,'mpmath')(R0,P0) for k,v in E.items()}
Gf=[[[sp.lambdify((r,p),Gam[k][i][j],'mpmath')(R0,P0) for j in range(4)] for i in range(4)] for k in range(4)]
gm=sp.lambdify((r,p),g,'mpmath')(R0,P0)
vf={k:[sp.lambdify((r,p),x,'mpmath') for x in u] for k,u in [('l',l),('n',n),('m',m)]}
c=mp.conj
def vec(k,R,PP):
    if k=='mb': return [c(f(R,PP)) for f in vf['m']]
    return [f(R,PP) for f in vf[k]]
W={'l':(1,1),'n':(-1,-1),'m':(1,-1),'mb':(-1,1)}
def covder(k,u):
    U=vec(u,R0,P0); out=[]
    for Aa in range(4):
        dd=U[2]*mp.diff(lambda x:vec(k,x,P0)[Aa],R0)+U[3]*mp.diff(lambda y:vec(k,R0,y)[Aa],P0)
        dd+=sum(Gf[Aa][Bb][Cc]*U[Bb]*vec(k,R0,P0)[Cc] for Bb in range(4) for Cc in range(4)); out.append(dd)
    return out
e,b,al,ga=sc['eps'],sc['beta'],sc['alp'],sc['gam']
ops={'Th':('l',lambda p_,q_:p_*e+q_*c(e)),'Thp':('n',lambda p_,q_:p_*ga+q_*c(ga)),'Ed':('m',lambda p_,q_:p_*b+q_*c(al)),'Edp':('mb',lambda p_,q_:p_*al+q_*c(b))}
dot=lambda x,y: sum(gm[i,j]*x[i]*y[j] for i in range(4) for j in range(4))
V0={k:vec(k,R0,P0) for k in W}
rp,tp=-sc['mu'],-sc['pi']; rho,tau=sc['rho'],sc['tau']
# predicted table (item 93(b))
pred={('Th','l'):{},('Th','n'):{'m':-tp,'mb':-c(tp)},('Th','m'):{'l':-c(tp)},('Th','mb'):{'l':-tp},
      ('Thp','l'):{'m':-c(tau),'mb':-tau},('Thp','n'):{},('Thp','m'):{'n':-tau},('Thp','mb'):{'n':-c(tau)},
      ('Ed','l'):{'m':-c(rho)},('Ed','n'):{'m':-rp},('Ed','m'):{},('Ed','mb'):{'l':-rp,'n':-c(rho)},
      ('Edp','l'):{'mb':-rho},('Edp','n'):{'mb':-c(rp)},('Edp','m'):{'l':-c(rp),'n':-rho},('Edp','mb'):{}}
worst=0
for op,(u,wt) in ops.items():
    for k in W:
        p_,q_=W[k]; Dv=covder(k,u); ww=wt(p_,q_); Dv=[Dv[i]-ww*V0[k][i] for i in range(4)]
        comps={'l':dot(Dv,V0['n']),'n':dot(Dv,V0['l']),'m':-dot(Dv,V0['mb']),'mb':-dot(Dv,V0['m'])}
        pr=pred[(op,k)]
        for kk in comps:
            err=abs(comps[kk]-pr.get(kk,0)); worst=max(worst,err)
print('tau, taubar-prime distinct here:',mp.nstr(tau,6),mp.nstr(-c(tp),6)); print('max deviation from the table:',mp.nstr(worst,5))

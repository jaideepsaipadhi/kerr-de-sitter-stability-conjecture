# Item 95: is the Killing-spinor vector xi of the (exact, nonlinear) Lambda-PD metric real up to a constant?
import sympy as sp, mpmath as mp, sys, time
t0=time.time(); mp.mp.dps=40
tau,sig,r,p=sp.symbols('tau sig r p',real=True); X=[tau,sig,r,p]
Lam=sp.Rational(1,20); a=sp.Rational(3,5)
def run(mv,nv,alv):
    m_,n_,al_=sp.Rational(mv),sp.Rational(nv),sp.Rational(alv)
    e_=1-Lam*a**2/3; k_=a**2
    P=k_+2*n_*p-e_*p**2+2*al_*m_*p**3-(al_**2*k_+Lam/3)*p**4
    Q=k_-2*m_*r+e_*r**2-2*al_*n_*r**3-(al_**2*k_+Lam/3)*r**4
    S=r**2+p**2; Om=1-al_*p*r
    w1=[1,-p**2,0,0]; w2=[1,r**2,0,0]
    g=sp.zeros(4)
    for i in range(4):
        for j in range(4):
            g[i,j]=(Q/S*w1[i]*w1[j]-P/S*w2[i]*w2[j])/Om**2
    g[2,2]=-S/Q/Om**2; g[3,3]=-S/P/Om**2
    gi=g.inv()
    V1=[r**2/S,-1/S,0,0]; V2=[p**2/S,1/S,0,0]
    al=Om*sp.sqrt(S/(2*Q))
    l=[al*V1[i]+al*(Q/S)*(1 if i==2 else 0) for i in range(4)]
    n=[al*V1[i]-al*(Q/S)*(1 if i==2 else 0) for i in range(4)]
    m=[Om/sp.sqrt(2)*(sp.sqrt(P/S)*(1 if i==3 else 0)+sp.I*sp.sqrt(S/P)*V2[i]) for i in range(4)]
    mb=[sp.conjugate(x) for x in m]
    Gam=[[[sum(gi[k,q]*(sp.diff(g[q,i],X[j])+sp.diff(g[q,j],X[i])-sp.diff(g[i,j],X[q]))/2 for q in range(4)) for j in range(4)] for i in range(4)] for k in range(4)]
    low=lambda u:[sum(g[i,j]*u[j] for j in range(4)) for i in range(4)]
    def cov(u):
        ul=low(u); return [[sp.diff(ul[i],X[j])-sum(Gam[k][i][j]*ul[k] for k in range(4)) for j in range(4)] for i in range(4)]
    Nl,Nn,Nmb=cov(l),cov(n),cov(mb)
    c3=lambda A,u,v: sum(A[i][j]*u[i]*v[j] for i in range(4) for j in range(4))
    E={'kap':c3(Nl,m,l),'rho':c3(Nl,m,mb),'tau':c3(Nl,m,n),'mu':-c3(Nn,mb,m),'pi':-c3(Nn,mb,l),
       'eps':(c3(Nl,n,l)+c3(Nmb,m,l))/2,'beta':(c3(Nl,n,m)+c3(Nmb,m,m))/2,'alp':(c3(Nl,n,mb)+c3(Nmb,m,mb))/2,'gam':(c3(Nl,n,n)+c3(Nmb,m,n))/2}
    F={k:sp.lambdify((r,p),v,'mpmath') for k,v in E.items()}
    vf={k:[sp.lambdify((r,p),x,'mpmath') for x in u] for k,u in [('l',l),('n',n),('m',m)]}
    R0,P0=mp.mpf(3),mp.mpf('0.36')
    sc={k:f(R0,P0) for k,f in F.items()}; c=mp.conj
    def D(Fn,u):
        U=[f(R0,P0) for f in vf[u]] if u!='mb' else [c(f(R0,P0)) for f in vf['m']]
        return U[2]*mp.diff(lambda x:Fn(x,P0),R0)+U[3]*mp.diff(lambda y:Fn(R0,y),P0)
    Thp=lambda Fn,pp,q: D(Fn,'n')-(pp*sc['gam']+q*c(sc['gam']))*Fn(R0,P0)
    Edp=lambda Fn,pp,q: D(Fn,'mb')-(pp*sc['alp']+q*c(sc['beta']))*Fn(R0,P0)
    rho_,tau_=sc['rho'],sc['tau']; rp_,tp_=-sc['mu'],-sc['pi']
    Fv=Thp(F['rho'],1,1)-Edp(F['tau'],1,-1)-(rho_*c(rp_)-tau_*c(tau_))
    Psi2=-Fv-2*Lam/6
    phi=Psi2**(-mp.mpf(1)/3)
    A,B,C,Dd=phi*rho_,-phi*rp_,phi*tau_,-phi*tp_
    Lv=[f(R0,P0) for f in vf['l']]; Nv=[f(R0,P0) for f in vf['n']]; Mv=[f(R0,P0) for f in vf['m']]; Mb=[c(z) for z in Mv]
    xi=[B*Lv[i]+A*Nv[i]-Dd*Mv[i]-C*Mb[i] for i in range(4)]
    print(f'm={mv} n={nv} alpha={alv}: |kappa|={mp.nstr(abs(sc["kap"]),3)} xi=',[mp.nstr(z,8) for z in xi],
          ' arg(xi^tau)-arg(xi^sig)=',mp.nstr((mp.arg(xi[0])-mp.arg(xi[1]))/mp.pi,10),'pi', round(time.time()-t0,1),flush=True)
for args in [(1,0,0),(1,'1/10',0),(1,0,'1/20'),(1,'1/10','1/20')]:
    run(*args)

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
Gf=[[[sp.lambdify((r,th),Gam[k][i][j],'mpmath') for j in range(4)] for i in range(4)] for k in range(4)]
G=[[[Gf[k][i][j](R0,T0) for j in range(4)] for i in range(4)] for k in range(4)]
gf=sp.lambdify((r,th),g,'mpmath'); gm=gf(R0,T0)
c=mp.conj
vf={k:[sp.lambdify((r,th),x,'mpmath') for x in u] for k,u in [('l',l),('n',n),('m',m)]}
def vecval(k,R,T):
    if k=='mb': return [c(f(R,T)) for f in vf['m']]
    return [f(R,T) for f in vf[k]]
W={'l':(1,1),'n':(-1,-1),'m':(1,-1),'mb':(-1,1)}
SC={k:f(R0,T0) for k,f in fn.items()} if False else None
sc={nm:sp.lambdify((r,th),e,'mpmath')(R0,T0) for nm,e in zip(['kap','rho','sig','tau','nu','pi','lam','mu','eps','beta','alp','gam'],[kap,rho,sig,tau,nu,pi,lam,mu,eps,beta,alp,gam])}
e,b,al_,ga=sc['eps'],sc['beta'],sc['alp'],sc['gam']
def covder(k,u):   # u^b nabla_b k^a at point
    U=vecval(u,R0,T0)
    out=[]
    for A in range(4):
        d=U[1]*mp.diff(lambda x:vecval(k,x,T0)[A],R0)+U[2]*mp.diff(lambda y:vecval(k,R0,y)[A],T0)
        d+=sum(G[A][B][C]*U[B]*vecval(k,R0,T0)[C] for B in range(4) for C in range(4))
        out.append(d)
    return out
ops={'Th':('l',lambda p,q:p*e+q*c(e)),'Thp':('n',lambda p,q:p*ga+q*c(ga)),'Ed':('m',lambda p,q:p*b+q*c(al_)),'Edp':('mb',lambda p,q:p*al_+q*c(b))}
dot=lambda x,y: sum(gm[i,j]*x[i]*y[j] for i in range(4) for j in range(4))
V0={k:vecval(k,R0,T0) for k in W}
rhop,taup=-sc['mu'],-sc['pi']
known={'rho':sc['rho'],'rhob':c(sc['rho']),'tau':sc['tau'],'taub':c(sc['tau']),'rp':rhop,'rpb':c(rhop),'tp':taup,'tpb':c(taup)}
for op,(u,wt) in ops.items():
    for k in W:
        p,q=W[k]
        D=covder(k,u); w=wt(p,q)
        D=[D[i]-w*V0[k][i] for i in range(4)]
        # components: X = (X.n) l + (X.l) n - (X.mb) m - (X.m) mb
        comps={'l':dot(D,V0['n']),'n':dot(D,V0['l']),'m':-dot(D,V0['mb']),'mb':-dot(D,V0['m'])}
        desc=[]
        for kk,vv in comps.items():
            if abs(vv)>1e-20:
                match=[nm for nm,val in known.items() if abs(vv-val)<1e-20]+['-'+nm for nm,val in known.items() if abs(vv+val)<1e-20]
                desc.append(kk+':'+(','.join(match) if match else mp.nstr(vv,6)))
        print(op,k,' '.join(desc),flush=True)

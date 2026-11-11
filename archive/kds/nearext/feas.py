import numpy as np
from scipy.integrate import solve_ivp
exec(open('/tmp/chk.py').read().split("for r0 in (0.9")[0])
def test(d,kap):
    t,A,B,z0,q1=d['t'],d['A'],d['B'],d['z0'],abs(d['q1'])
    mu=lambda z: z*(z-1)*(t-z); mup=lambda z: (z-1)*(t-z)+z*(t-z)-z*(z-1)
    q=lambda z: q1*(z-z0)+kap*(z-1)*(t-z); qp=lambda z: q1+kap*(t+1-2*z); qpp=-2*kap
    L=lambda z: A/z+B-3*z+(kap/z)*(2*q1*(z-z0)+kap*(z-1)*(t-z))
    # zero of q
    from scipy.optimize import brentq
    zt=brentq(q,1+1e-12,t-1e-9)
    zs=np.linspace(zt,t,4000)[1:-1]
    Lp=lambda z,h=1e-6*t: (L(z+h)-L(z-h))/(2*h)
    Cp=mup(zs)*q(zs)-2*mu(zs)*qp(zs); Dp=q(zs)*Lp(zs)-kap*mup(zs)
    cos_ok=(Cp.max()<=0) and (Dp.max()<=0) and L(t)>0
    # event side: D f=0 ODE from z=1
    def rhs(z,y):
        f,f1=y; muz=mu(z); a2=0.5*muz*qp(z); a1=0.5*(mup(z)*qp(z)+muz*qpp)+0.5*muz*qpp+q(z)*L(z)
        a0=0.5*(mup(z)*qpp)+q(z)*Lp(z)
        return [f1,-(a1*f1+a0*f)/a2]
    # Frobenius exponent at 1: approximate p from leading balance
    p=2*abs(q(1))*L(1)/(qp(1)*mup(1))
    e=1e-7*(zt-1); y0=[e**p,p*e**(p-1)]
    sol=solve_ivp(rhs,[1+e,zt],y0,rtol=1e-10,atol=1e-300,dense_output=True)
    zz=np.linspace(1+e,zt,3000); f,f1=sol.sol(zz)
    C=mup(zz)*f*q(zz)-2*mu(zz)*f*qp(zz)-mu(zz)*f1*q(zz)
    ev_ok=(f.min()>0) and (C.max()<=1e-9*abs(C).max()) and f1[-1]>=0
    return cos_ok, ev_ok, L(t), zt, p
for r0,frac in []:
    d=data(r0,1.0,10.0,2,frac)
    for kap in [0, -d['Lt']*d['t']/(2*d['q1']*(d['t']-d['z0']))*1.5, -d['Lt']*d['t']/(2*d['q1']*(d['t']-d['z0']))*5]:
        print(r0,frac,round(kap,4),test(d,kap))
from scipy.integrate import quad
def event_R3(d,kap,n=6000):
    t,A,B,z0,q1=d['t'],d['A'],d['B'],d['z0'],abs(d['q1'])
    mu=lambda z: z*(z-1)*(t-z); mup=lambda z: (z-1)*(t-z)+z*(t-z)-z*(z-1)
    q=lambda z: q1*(z-z0)+kap*(z-1)*(t-z); qp=lambda z: q1+kap*(t+1-2*z); qpp=-2*kap
    L=lambda z: A/z+B-3*z+(kap/z)*(2*q1*(z-z0)+kap*(z-1)*(t-z))
    Lp=lambda z: -A/z**2-3+kap*(2*q1*z0+... if False else 0) 
    h=lambda z,e=1e-7: (L(z+e)-L(z-e))/(2*e)
    from scipy.optimize import brentq
    zt=brentq(q,1+1e-12,t-1e-9)
    hC=lambda z: (2*mu(z)*qp(z)-mup(z)*q(z))/(-mu(z)*q(z))
    zz=np.linspace(1+1e-9,zt-1e-9*(zt-1),n)
    # log F_C via cumulative integral (F_C ~ (z-1) near 1): integrate hC - 1/(z-1) then add log(z-1)
    g=np.array([hC(z)-1/(z-1) for z in zz]); lF=np.log(zz-1)+np.concatenate([[0],np.cumsum((g[1:]+g[:-1])/2*np.diff(zz))])
    F=np.exp(lF-lF.max()); F1=F*np.array([hC(z) for z in zz])
    # D(F)= q(F L)' + 1/2 (mu (F q'))''  numerically
    G=mu(zz)*np.gradient(F*qp(zz),zz)
    D=q(zz)*np.gradient(F*L(zz),zz)+0.5*np.gradient(G,zz)
    sgn=np.sign(D[5:-5]); ch=np.where(np.diff(sgn)!=0)[0]
    if len(ch)==0: return ('no sign change', D[10]<0)
    i=ch[0]+5; zh=zz[i]
    def rhs(z,y):
        f,f1=y; muz=mu(z)
        a2=0.5*muz*qp(z); a1=0.5*(mup(z)*qp(z)+muz*qpp)+0.5*muz*qpp+q(z)*L(z); a0=0.5*mup(z)*qpp+q(z)*h(z)
        return [f1,-(a1*f1+a0*f)/a2]
    sol=solve_ivp(rhs,[zh,zt],[F[i],F1[i]],rtol=1e-10,atol=1e-300)
    f,f1=sol.y
    return (len(ch), round(zh,4), bool(f.min()>0), bool(f1[-1]>=0))
for r0,frac in []:
    d=data(r0,1.0,10.0,2,frac)
    k0=-d['Lt']*d['t']/(2*d['q1']*(d['t']-d['z0']))
    for kap in [0,1.5*k0 if k0>0 else 0,3*k0 if k0>0 else 0]:
        print(r0,frac,round(kap,3),event_R3(d,kap))
for r0,frac in [(0.3,0.5),(0.5,0.5)]:
    d=data(r0,1.0,10.0,2,frac); print('sanity',r0,d['B']>d['floor'],event_R3(d,0))
    t,A,B,z0=d['t'],d['A'],d['B'],d['z0']
    x0=z0-1; p=2*x0*(A+B-3)/(t-1); print(' p=',p,' R2cond',2*x0**2*(A/z0+3)<t-1)

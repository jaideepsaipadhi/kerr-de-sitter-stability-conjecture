import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
exec(open('/tmp/chk.py').read().split("for r0 in (0.9")[0])
def event(d,kap,n=20000):
    t,A,B,z0,q1=d['t'],d['A'],d['B'],d['z0'],abs(d['q1'])
    mu=lambda z: z*(z-1)*(t-z); mup=lambda z: (z-1)*(t-z)+z*(t-z)-z*(z-1); mupp=lambda z: 2*(t+1)-6*z
    q=lambda z: q1*(z-z0)+kap*(z-1)*(t-z); qp=lambda z: q1+kap*(t+1-2*z); qpp=-2*kap
    L=lambda z: A/z+B-3*z+(kap/z)*(2*q1*(z-z0)+kap*(z-1)*(t-z))
    e=1e-7
    Lp=lambda z: (L(z+e*z)-L(z-e*z))/(2*e*z)
    zt=brentq(q,1+1e-12,t-1e-9)
    hC=lambda z: (2*mu(z)*qp(z)-mup(z)*q(z))/(-mu(z)*q(z))
    hCp=lambda z: (hC(z+e*(zt-1))-hC(z-e*(zt-1)))/(2*e*(zt-1))
    def Dn(z,h,hp):  # D^f / f
        mq=mup(z)*qp(z)+mu(z)*qpp
        return q(z)*(L(z)*h+Lp(z))+0.5*(mq*h+mu(z)*qp(z)*(hp+h*h)+mup(z)*qpp+mu(z)*qpp*h)
    zz=1+(zt-1)*np.linspace(1e-4,1-1e-4,n)
    dv=np.array([Dn(z,hC(z),hCp(z)) for z in zz])
    neg0=dv[0]<0
    ch=np.where(np.diff(np.sign(dv))!=0)[0]
    if len(ch)==0: return ('D(F_C) one sign', 'neg' if neg0 else 'pos')
    i=ch[0]; zh=zz[i]
    def rhs(z,y):
        f,f1=y; a2=0.5*mu(z)*qp(z); mq=mup(z)*qp(z)+mu(z)*qpp
        a1=0.5*mq+0.5*mu(z)*qpp+q(z)*L(z); a0=0.5*mup(z)*qpp+q(z)*Lp(z)
        return [f1,-(a1*f1+a0*f)/a2]
    sol=solve_ivp(rhs,[zh,zt-1e-9*(zt-1)],[1.0,hC(zh)],rtol=1e-10,atol=1e-14,dense_output=True)
    z2=np.linspace(zh,sol.t[-1],2000); f,f1=sol.sol(z2)
    Cok=bool(np.all(f1/f<=np.array([hC(z) for z in z2])*(1+1e-8)+1e-12))
    return dict(changes=len(ch),neg_start=bool(neg0),zh=round(zh,4),fpos=bool(f.min()>0),seam=bool(f1[-1]>=0),C=Cok)
if __name__=='__main__':
  for r0,frac in [(0.3,0.5),(0.99,0.99),(0.99,0.9),(0.9,0.99),(0.9,0.9),(0.99,0.5)]:
    d=data(r0,1.0,10.0,2,frac)
    k0=-d['Lt']*d['t']/(2*d['q1']*(d['t']-d['z0']))
    ks=[0] if k0<=0 else [0,1.5*k0,3*k0]
    for kap in ks: print(r0,frac,round(kap,3),event(d,kap))

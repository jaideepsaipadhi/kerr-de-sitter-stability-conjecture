import numpy as np
from scipy.integrate import solve_ivp
rng=np.random.default_rng(31)
stats={'n':0,'single_change':0,'d_neg_near1':0,'ok':0,'failC':0,'failpos':0,'failseam':0,'multi':0}
ex=[]
for _ in range(400):
    t=1+10**rng.uniform(-2,4); zf=rng.uniform(0.02,0.98); z0=1+zf*(t-1)
    B=(5.5*t-1.25)*10**rng.uniform(0,2); A=B*10**rng.uniform(-2,5)
    if 2*(z0-1)*(A+B-3)/(t-1)<=1 or 2*(z0-1)**2*(A/z0+3)<t-1: continue   # R3 only
    mu=lambda z: z*(z-1)*(t-z); mup=lambda z: -3*z*z+2*(1+t)*z-t; mupp=lambda z: -6*z+2*(1+t)
    L=lambda z: A/z+B-3*z; Lp=lambda z: -A/z**2-3
    H=lambda z: 2*mu(z)+(z0-z)*mup(z); Hp=lambda z: 2*mup(z)-mup(z)+(z0-z)*mupp(z)
    hC=lambda z: H(z)/((z0-z)*mu(z))
    def hCp(z):
        s_=z0-z; Hp_=mup(z)+s_*mupp(z); return (Hp_*s_*mu(z)-H(z)*(-mu(z)+s_*mup(z)))/(s_*mu(z))**2
    d=lambda z: 0.5*(mup(z)*hC(z)+mu(z)*(hCp(z)+hC(z)**2))-(z0-z)*(L(z)*hC(z)+Lp(z))
    zs=1+(z0-1)*np.linspace(1e-4,1-1e-4,4000); dv=np.array([d(z) for z in zs])
    stats['n']+=1
    sc=np.where(np.diff(np.sign(dv))!=0)[0]
    if dv[0]<0: stats['d_neg_near1']+=1
    if len(sc)!=1 or dv[0]>=0 or dv[-1]<=0: stats['multi']+=1; ex.append((t,zf,A/B,B/(5.5*t),len(sc),dv[0]>0)); continue
    stats['single_change']+=1
    a_,b_=zs[sc[0]],zs[sc[0]+1]
    for _i in range(80):
        m_=(a_+b_)/2
        if d(m_)<0: a_=m_
        else: b_=m_
    zh=b_
    # phase 2: (mu f')' = 2 s (L f)'  ->  mu f'' + (mu' - 2 s L) f' - 2 s L' f = 0
    def rhs(z,Y):
        f,fp=Y; s=z0-z
        return [fp, (-(mup(z)-2*s*L(z))*fp+2*s*Lp(z)*f)/mu(z)]
    sol=solve_ivp(rhs,[zh,z0*(1-1e-9)+1e-9*zh if False else z0-1e-9*(z0-1)],[1.0,hC(zh)],rtol=1e-10,atol=1e-13,dense_output=True)
    zz=np.linspace(zh,sol.t[-1],2000); F=sol.sol(zz)
    if np.any(F[0]<=0): stats['failpos']+=1; continue
    hh=F[1]/F[0]; hc=np.array([hC(z) for z in zz])
    viol=np.max((hh-hc)/np.abs(hc)); stats.setdefault("maxviol",-9); stats["maxviol"]=max(stats["maxviol"],viol)
    if viol>1e-3: stats["failC"]+=1; ex.append((t,zf,A/B,viol,np.argmax((hh-hc)/np.abs(hc))/len(hh))); continue
    if F[1][-1]<=0: stats['failseam']+=1; continue
    stats['ok']+=1
print(stats); print(ex[:8])

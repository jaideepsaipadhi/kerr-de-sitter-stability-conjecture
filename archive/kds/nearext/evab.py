import numpy as np, sympy as sp, time
z,t_,z0_,A_,B_,S=sp.symbols('z t z0 A B S')
mu=z*(z-1)*(t_-z); s=z0_-z
F=mu/s**2; L=A_/z+B_-3*z
Dt=sp.diff(mu*sp.diff(F,z),z)-2*s*sp.diff(L*F,z)
Nexpr=sp.simplify(Dt*s**4)
Nf=sp.lambdify((z,t_,z0_,A_,B_),Nexpr,'numpy')
rng=np.random.default_rng(5); bad=[];N=0;t0=time.time()
from scipy.integrate import solve_ivp
while time.time()-t0<40:
    t=10**rng.uniform(0.01,3); z0=1+(t-1)*rng.uniform(0.001,0.999); A=10**rng.uniform(-3,6)
    Bf=5.5*t-1.25; B=Bf-10**rng.uniform(-2,np.log10(Bf+1e4))
    Lmin=min(A/zz+B-3*zz for zz in np.linspace(1,z0,50))
    if A/z0+B-3*z0<=0: continue
    x0=z0-1; p=2*x0*(A+B-3)/(t-1)
    if p<=1: continue  # R1 handled (needs L>0? check below)
    zz=1+x0*np.linspace(1e-6,1-1e-6,4000); Nv=Nf(zz,t,z0,A,B)
    ch=np.sum(np.diff(np.sign(Nv))!=0)
    N+=1
    if ch!=1 or Nv[0]>=0: bad.append((round(t,2),round((z0-1)/(t-1),3),round(A,2),round(B,2),round(Lmin,2),int(ch))); continue
    # seam: K(zhat)>0 test
    i=np.where(np.diff(np.sign(Nv))!=0)[0][0]; zh=zz[i]; sh=z0-zh
    muh=zh*(zh-1)*(t-zh); mup=(zh-1)*(t-zh)+zh*(t-zh)-zh*(zh-1); H=2*muh+sh*mup
    K=H-sh**3*(A/zh**2+3)
    if K<=0: bad.append(('seamK',round(t,2),round((z0-1)/(t-1),3),round(A,2),round(B,2),round(Lmin,2)))
print('N',N,'bad',len(bad)); [print(b) for b in bad[:12]]

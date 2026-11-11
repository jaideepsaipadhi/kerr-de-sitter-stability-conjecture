import numpy as np
from scipy.integrate import solve_ivp
exec(open('/tmp/evab.py').read().split("rng=")[0])
def seam(t,z0,A,B):
    x0=z0-1; zz=1+x0*np.linspace(1e-6,1-1e-6,4000); Nv=Nf(zz,t,z0,A,B)
    i=np.where(np.diff(np.sign(Nv))!=0)[0][0]; zh=zz[i]
    mu=lambda z: z*(z-1)*(t-z); mup=lambda z:(z-1)*(t-z)+z*(t-z)-z*(z-1)
    L=lambda z:A/z+B-3*z; Lp=lambda z:-A/z**2-3
    def rhs(z,y):
        f,f1=y; s=z0-z; return [f1,(-(mup(z)-2*s*L(z))*f1+2*s*Lp(z)*f)/mu(z)]
    sh=z0-zh; f1=(2*mu(zh)+sh*mup(zh))/sh**3; f=mu(zh)/sh**2
    sol=solve_ivp(rhs,[zh,z0],[1.0,f1/f],rtol=1e-11,atol=1e-14)
    return sol.y[0].min()>0, sol.y[1][-1]
rng=np.random.default_rng(5); bad=0; n=0; worst=[]
import time; t0=time.time()
while time.time()-t0<45:
    t=10**rng.uniform(0.01,3); z0=1+(t-1)*rng.uniform(0.001,0.999); A=10**rng.uniform(-3,6)
    Bf=5.5*t-1.25; B=Bf-10**rng.uniform(-2,np.log10(Bf+1e4))
    if A/z0+B-3*z0<=0: continue
    x0=z0-1; p=2*x0*(A+B-3)/(t-1)
    if p<=1: continue
    zz=1+x0*np.linspace(1e-6,1-1e-6,2000); Nv=Nf(zz,t,z0,A,B)
    i=np.where(np.diff(np.sign(Nv))!=0)[0]
    zh=zz[i[0]]; sh=z0-zh; H=2*zh*(zh-1)*(t-zh)+sh*((zh-1)*(t-zh)+zh*(t-zh)-zh*(zh-1)); K=H-sh**3*(A/zh**2+3)
    if K>0: continue
    n+=1; ok,fs=seam(t,z0,A,B)
    if not(ok and fs>0): bad+=1; worst.append((round(t,2),round(x0/(t-1),3),round(A,1),round(B,1),fs))
print('K-fail cases',n,'true seam fails',bad); print(worst[:8])

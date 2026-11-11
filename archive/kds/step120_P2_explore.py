import numpy as np
from scipy.integrate import solve_ivp
rng=np.random.default_rng(51); rows=[]
for _ in range(600):
    t=1+10**rng.uniform(-2,4); zf=rng.uniform(0.02,0.98); z0=1+zf*(t-1)
    B=(5.5*t-1.25)*10**rng.uniform(0,2); A=B*10**rng.uniform(-2,5)
    if 2*(z0-1)*(A+B-3)/(t-1)<=1 or 2*(z0-1)**2*(A/z0+3)<t-1: continue
    mu=lambda z: z*(z-1)*(t-z); mup=lambda z: -3*z*z+2*(1+t)*z-t; mupp=lambda z: -6*z+2*(1+t)
    L=lambda z: A/z+B-3*z; Lp=lambda z: -A/z**2-3; Lpp=lambda z: 2*A/z**3
    H=lambda z: 2*mu(z)+(z0-z)*mup(z); hC=lambda z: H(z)/((z0-z)*mu(z))
    def hCp(z):
        s_=z0-z; return ((mup(z)+s_*mupp(z))*s_*mu(z)-H(z)*(-mu(z)+s_*mup(z)))/(s_*mu(z))**2
    d=lambda z: 0.5*(mup(z)*hC(z)+mu(z)*(hCp(z)+hC(z)**2))-(z0-z)*(L(z)*hC(z)+Lp(z))
    zs=1+(z0-1)*np.linspace(1e-4,1-1e-4,3000); dv=np.array([d(z) for z in zs]); i=np.where(np.diff(np.sign(dv))!=0)[0][0]
    a_,b_=zs[i],zs[i+1]
    for _k in range(80):
        m_=(a_+b_)/2
        if d(m_)<0: a_=m_
        else: b_=m_
    zh=b_
    def rhs(z,Y):  # Riccati for h
        h=Y[0]; s=z0-z; return [ (2/mu(z))*(s*(L(z)*h+Lp(z))-0.5*mup(z)*h)-h*h ]
    sol=solve_ivp(rhs,[zh,z0],[hC(zh)],rtol=1e-11,atol=1e-14,dense_output=True)
    zz=np.linspace(zh,z0,1500); hh=sol.sol(zz)[0]
    g=np.array([-Lp(z)/L(z) for z in zz])
    rows.append((t,zf,A/B,hh[-1]/g[-1],np.min(hh/g),(zh-1)/(z0-1)))
R=np.array(rows)
print('n',len(R),' min h(z0)/g(z0)',R[:,3].min(),' min over phase2 of h/g',R[:,4].min(),' zhat fraction range',R[:,5].min(),R[:,5].max())
i=np.argmin(R[:,4]); print('worst',R[i])
rng=np.random.default_rng(52); res={}
cands={'0.5hC':lambda hc,g: 0.5*hc,'0.3hC':lambda hc,g:0.3*hc,'0.7hC':lambda hc,g:0.7*hc,'2g':lambda hc,g:2*g,'3g':lambda hc,g:3*g,'g':lambda hc,g:g,
       'hC-2/s':None}
cnt={k:0 for k in cands}; n=0
for _ in range(600):
    t=1+10**rng.uniform(-2,4); zf=rng.uniform(0.02,0.98); z0=1+zf*(t-1)
    B=(5.5*t-1.25)*10**rng.uniform(0,2); A=B*10**rng.uniform(-2,5)
    if 2*(z0-1)*(A+B-3)/(t-1)<=1 or 2*(z0-1)**2*(A/z0+3)<t-1: continue
    mu=lambda z: z*(z-1)*(t-z); mup=lambda z: -3*z*z+2*(1+t)*z-t; mupp=lambda z: -6*z+2*(1+t)
    L=lambda z: A/z+B-3*z; Lp=lambda z: -A/z**2-3; Lpp=lambda z: 2*A/z**3
    H=lambda z: 2*mu(z)+(z0-z)*mup(z); hC=lambda z: H(z)/((z0-z)*mu(z))
    def hCp(z):
        s_=z0-z; return ((mup(z)+s_*mupp(z))*s_*mu(z)-H(z)*(-mu(z)+s_*mup(z)))/(s_*mu(z))**2
    d=lambda z: 0.5*(mup(z)*hC(z)+mu(z)*(hCp(z)+hC(z)**2))-(z0-z)*(L(z)*hC(z)+Lp(z))
    zs=1+(z0-1)*np.linspace(1e-4,1-1e-4,3000); dv=np.array([d(z) for z in zs]); i=np.where(np.diff(np.sign(dv))!=0)[0][0]; zh=zs[i+1]
    zz=np.linspace(zh,z0-1e-9*(z0-1),800); n+=1
    hc=np.array([hC(z) for z in zz]); hcp=np.array([hCp(z) for z in zz]); g=np.array([-Lp(z)/L(z) for z in zz])
    gp=g*g-np.array([Lpp(z)/L(z) for z in zz])
    s=z0-zz; MU=np.array([mu(z) for z in zz]); MUP=np.array([mup(z) for z in zz]); LL=np.array([L(z) for z in zz]); LP=np.array([Lp(z) for z in zz])
    for k,f in cands.items():
        if f is None: continue
        kk=f(hc,g); kp=f(hcp,gp)
        D=0.5*(MUP*kk+MU*(kp+kk*kk))-s*(LL*kk+LP)
        if np.all(D<=0): cnt[k]+=1
print('n',n,cnt)

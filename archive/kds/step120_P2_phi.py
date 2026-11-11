import numpy as np
rng=np.random.default_rng(61); n=0; bad=0; phiok=0
for _ in range(40000):
    t=1+10**rng.uniform(-2,5); zf=rng.uniform(0,1); z0=1+zf*(t-1); x0=z0-1
    B=(5.5*t-1.25)*10**rng.uniform(0,3); A=B*10**rng.uniform(-3,6)
    if 2*x0*(A+B-3)/(t-1)<=1 or 2*x0**2*(A/z0+3)<t-1: continue
    z=1+x0*np.linspace(1e-6,1-1e-6,4000); s=z0-z
    mu=z*(z-1)*(t-z); mup=-3*z*z+2*(1+t)*z-t; L=A/z+B-3*z
    H=2*mu+s*mup; J=H-2*s*s*L
    S=s; mu0=z0*x0*(t-z0); mup0=-3*z0**2+2*(1+t)*z0-t
    c=[6*mu0**2,-6*mu0*mup0,4*(z0-1)*(z0-t)*(A+B*z0)+6*t**2*z0**2-6*t**2*z0+t**2-8*t*z0**3+12*t*z0**2-6*t*z0+3*z0**4-8*z0**3+6*z0**2,
       2*(A*(t+1-2*z0)+B*mup0+3*z0*(-3*t*z0+2*t+4*z0**2-3*z0)),0,2*(B+4*t-15*z0+4),15]
    N=sum(c[k]*S**k for k in range(7)); n+=1
    if np.any((J<=0)&(N>=0)): bad+=1
    # direct: J at zhat
    i=np.where(np.diff(np.sign(N))!=0)[0]
    if len(i): phiok+= J[i[0]]>0 or J[i[0]+1]>0
print('n',n,'violations of [J<=0 => N<0]:',bad,' J(zhat)>0 count',phiok)

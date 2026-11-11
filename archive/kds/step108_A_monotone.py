import numpy as np
rng=np.random.default_rng(3); bad=0; tot=0; inc=dec=0
def F(t,z0,A,z):
    B=5.5*t-1.25; L=A/z+B-3*z; N=A/z**2+3; x=z-1; x0=z0-1; s=z0-z
    mu=z*x*(t-z); mup=-3*z*z+2*(1+t)*z-t; c=np.where(mup>0,0.85,17/22)
    L0=A/z0+B-3*z0; w=8*L0*s*s/(z0*x0*(t-z0))
    ok=(w>=10)&(2*x0*(A+B-3)/(t-1)>1)&(2*x0**2*(A/z0+3)>=t-1)
    return c*mup+2.295*mu*N/L-0.5455*s*L-0.225*(t-z)*L/(x0*N), ok
for _ in range(20000):
    t=1+10**rng.uniform(-2,5); z0=1+rng.uniform(0,1)*(t-1); z=1+rng.uniform(0,1)*(z0-1)
    As=10**np.linspace(-4,8,200)*(5.5*t)
    f,ok=F(t,z0,As,z)
    if ok.sum()<3: continue
    tot+=1; fo=f[ok]; d=np.diff(fo)
    if np.all(d<=0): dec+=1
    elif np.all(d>=0): inc+=1
    else:
        # unimodal check: max at an endpoint?
        if np.argmax(fo) in (0,len(fo)-1): bad+=0
        else: bad+=1
print('tot',tot,'monotone dec',dec,'inc',inc,'interior max',bad)

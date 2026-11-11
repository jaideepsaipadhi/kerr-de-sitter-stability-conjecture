import numpy as np
rng=np.random.default_rng(11)
def F(t,z0,z,A,B):
    x=z-1;x0=z0-1;s=z0-z;mu=z*x*(t-z);mup=-3*z*z+2*(1+t)*z-t
    L=A/z+B-3*z; g=(A/z**2+3)/L; v=1-x/x0; q=0.45*v/(x*z); qp=0.45*(-(1/x0)/(x*z)-v*(2*z-1)/(x*z)**2)
    h=np.where(mup>0,1.7*g+q,17/11*g+q)
    val=0.5*mup*h+0.85*mu*g*g+0.5*mu*(1.7*g+q)**2+0.5*mu*qp-0.5455*s*L*g-s*L*q
    L0=A/z0+B-3*z0; w=8*L0*s*s/(z0*x0*(t-z0)); ok=(w>=10)&(2*x0*(A+B-3)/(t-1)>1)&(2*x0**2*(A/z0+3)>=t-1)
    return val,ok
stats={'dec':0,'inc':0,'mixed':0,'mixed_interior_max':0}
for _ in range(20000):
    t=1+10**rng.uniform(-2,5); z0=1+rng.uniform(0,1)*(t-1); z=1+rng.uniform(0,1)*(z0-1)
    A=(5.5*t)*10**rng.uniform(-2,6); Bs=(5.5*t-1.25)*10**np.linspace(0,3,150)
    f,ok=F(t,z0,z,A,Bs)
    if ok.sum()<3: continue
    fo=f[ok]; d=np.diff(fo)
    if np.all(d<=0): stats['dec']+=1
    elif np.all(d>=0): stats['inc']+=1
    else:
        stats['mixed']+=1
        if 0<np.argmax(fo)<len(fo)-1: stats['mixed_interior_max']+=1
print(stats)
# slice 1: B = floor, vary A (i.e. L): where is the max?
st={'Lmin':0,'interior':0,'n':0}
st2={'Bfloor':0,'interior':0,'top':0,'n':0}
for _ in range(20000):
    t=1+10**rng.uniform(-2,5); z0=1+rng.uniform(0,1)*(t-1); z=1+rng.uniform(0,1)*(z0-1)
    B=5.5*t-1.25; As=(5.5*t)*np.concatenate([[0],10**np.linspace(-4,7,300)])
    f,ok=F(t,z0,z,As,B)
    if ok.sum()>=3:
        st['n']+=1; i=np.argmax(np.where(ok,f,-np.inf)); first=np.argmax(ok)
        st['Lmin' if i==first else 'interior']+=1
    # slice 2: A = A_min(B) (smallest feasible A for each B), vary B
    Bs=(5.5*t-1.25)*10**np.linspace(0,3,150); vals=[]
    x0=z0-1; s=z0-z; mu0=z0*x0*(t-z0)
    for Bv in Bs:
        Am=max(0, z0*(10*mu0/(8*s*s)-Bv+3*z0), (t-1)/(2*x0)-Bv+3+1e-9, z0*((t-1)/(2*x0*x0)-3))
        fv,okv=F(t,z0,z,np.array([Am*(1+1e-12)+1e-12]),Bv); vals.append(fv[0] if okv[0] else -np.inf)
    vals=np.array(vals)
    if np.isfinite(vals).sum()>=3:
        st2['n']+=1; i=np.argmax(vals); fi=np.argmax(np.isfinite(vals))
        st2['Bfloor' if i==fi else ('top' if i==len(vals)-1 else 'interior')]+=1
print(st); print(st2)

import numpy as np
rng=np.random.default_rng(22); bad=0; n=0; worst=9
for _ in range(40000):
    t=1+10**rng.uniform(-2,5); z0=1+rng.uniform(0,1)*(t-1); z=1+rng.uniform(0,1)*(z0-1)
    x=z-1;x0=z0-1;s=z0-z;mu=z*x*(t-z);mup=-3*z*z+2*(1+t)*z-t;mu0=z0*x0*(t-z0)
    v=1-x/x0;q=0.45*v/(x*z);qp=0.45*(-(1/x0)/(x*z)-v*(2*z-1)/(x*z)**2)
    c=1.7 if mup>0 else 17/11
    L0s=1.25*mu0/s**2; kap=s/z; L1=(z0/z)*(L0s+3*z0)-3*z
    g0=(1+1/kap)/z; C=(L1/kap-6*z)/z
    G1c=0.5*mup*c+1.7*mu*q; G0c=0.5*mup*q+0.5*mu*q*q+0.5*mu*qp
    a2=2.295*mu*C*C; a1=-C*(4.59*mu*g0+G1c); E=s*(0.5455*g0+q)
    fl=5.5*t-1.25; Bhi=L0s+3*z0
    Bpp=(z0*(L0s+3*z0)-3-(t-1)/(2*x0))/(z0-1) ; Bc2=L0s+3*z0-((t-1)/(2*x0*x0)-3)
    Bmax=min(Bhi,Bpp,Bc2)
    if Bmax<=fl: continue
    Llo=L1-kap*Bmax; Lhi=L1-kap*fl
    if Llo<=0: continue
    us=np.linspace(1/Lhi,1/Llo,2000); n+=1
    fp=2*a2*us+a1+E/us**2; fpp=2*a2-2*E/us**3
    m=fpp<=0
    if m.any():
        mm=fp[m].min()/ (abs(a1)+2*a2*us[m].max()+E/us[m].min()**2)
        worst=min(worst,mm)
        if mm<=0: bad+=1
print('lines',n,'violations of [f\'>0 where f\'\'<=0]',bad,'worst normalised f\' there',worst)
rng=np.random.default_rng(23); bad=0;n=0;worst=9
for _ in range(100000):
    t=1+10**rng.uniform(-2,5); z0=1+rng.uniform(0,1)*(t-1); z=1+rng.uniform(0,1)*(z0-1)
    x=z-1;x0=z0-1;s=z0-z;mu=z*x*(t-z);mup=-3*z*z+2*(1+t)*z-t;mu0=z0*x0*(t-z0)
    v=1-x/x0;q=0.45*v/(x*z);qp=0.45*(-(1/x0)/(x*z)-v*(2*z-1)/(x*z)**2)
    c=1.7 if mup>0 else 17/11
    L0s=1.25*mu0/s**2; kap=s/z; L1=(z0/z)*(L0s+3*z0)-3*z
    g0=(1+1/kap)/z; C=(L1/kap-6*z)/z
    G1c=0.5*mup*c+1.7*mu*q
    a2=2.295*mu*C*C; a1=-C*(4.59*mu*g0+G1c); E=s*(0.5455*g0+q)
    fl=5.5*t-1.25; Bhi=L0s+3*z0
    Bpp=(z0*(L0s+3*z0)-3-(t-1)/(2*x0))/(z0-1); Bc2=L0s+3*z0-((t-1)/(2*x0*x0)-3)
    if min(Bhi,Bpp,Bc2)<=fl: continue
    Lhi=L1-kap*fl; ulo=1/Lhi
    if E<=a2*ulo**3: continue
    n+=1
    if a1>=0: continue
    val=(27*a2*a2*E+a1**3)/(27*a2*a2*E+abs(a1)**3)
    worst=min(worst,val); bad+= val<=0
print('C*: lines',n,'violations',bad,'worst',worst)

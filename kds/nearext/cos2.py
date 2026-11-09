import numpy as np
rng=np.random.default_rng(0)
def cos_ok(t,z0,A,B,kap,n=4000):
    z=z0+(t-z0)*np.linspace(1e-6,1,n)
    mu=z*(z-1)*(t-z); mup=(z-1)*(t-z)+z*(t-z)-z*(z-1)
    q1=1.0  # scale: conditions homogeneous? q1 enters L-tilde -> keep q1 symbolic via kap/q1 ratio
    qt=q1*(z-z0)+kap*(z-z0)*(t-z); qtp=q1+kap*(t+z0-2*z)
    Lt=A/z+B-3*z+kap*(z-z0)*(2*q1*(z-z0)+kap*(z-z0)*(t-z))/(z*(z-1))
    Ltp=np.gradient(Lt,z)
    C=mup*qt-2*mu*qtp; D=qt*Ltp-kap*mup
    return C.max()<=1e-9*abs(C).max(), D[1:-1].max()<=1e-9*abs(D).max(), Lt[-1]>0

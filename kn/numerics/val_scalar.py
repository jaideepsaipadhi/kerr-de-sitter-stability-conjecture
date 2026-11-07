import numpy as np
from scipy.integrate import solve_ivp
s=-2; lam=4.0; M=1
def F(r,y,w):
    R=y[0]+1j*y[1]; dR=y[2]+1j*y[3]; D=r*r-2*r; K=r*r*w
    d2=-((s+1)*(2*r-2)/D*dR+((K*K-2j*s*(r-1)*K)/D+4j*s*w*r-lam)/D*R)
    return [dR.real,dR.imag,d2.real,d2.imag]
def W(w,rbig=60.,rmid=3.,x0=1e-5):
    q=-s-1j*w*4  # (r-2)^{-s-i w (r+^2)/(r+-r-)} = -s - 4iw
    r0=2+x0; y=[x0**q,0,q*x0**(q-1),0]; y=[np.real(y[0]),np.imag(y[0]),np.real(y[2]),np.imag(y[2])]
    a=solve_ivp(F,(r0,rmid),y,args=(w,),rtol=1e-10,atol=1e-30).y[:,-1]
    rs=rbig+2*np.log(rbig/2-1); ph=np.exp(1j*w*rs)*rbig**3; dph=ph*(3/rbig+1j*w*rbig/(rbig-2))
    b=solve_ivp(F,(rbig,rmid),[ph.real,ph.imag,dph.real,dph.imag],args=(w,),rtol=1e-10,atol=1e-30).y[:,-1]
    A=np.array([[a[0]+1j*a[1],b[0]+1j*b[1]],[a[2]+1j*a[3],b[2]+1j*b[3]]]); A/=np.linalg.norm(A,axis=0)
    return abs(np.linalg.det(A))
for w in [0.37367168-0.08896232j,0.40-0.089j,0.35-0.089j,0.3737-0.06j,0.3737-0.12j]: print(w,W(w))

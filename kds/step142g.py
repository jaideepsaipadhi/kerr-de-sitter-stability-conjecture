exec(open('step141b.py').read().split("for c in (0.5")[0])
from scipy.integrate import cumulative_trapezoid as ct
ok=0
for P in out:
    tt,zz,AA,BB,Sh=P; zh=zz-Sh; mu_=zh*(zh-1)*(tt-zh); mup=(zh-1)*(tt-zh)+zh*(tt-zh)-zh*(zh-1); w0=(2*mu_+Sh*mup)/Sh
    q=np.sort(zh+Sh*np.concatenate([[0],np.geomspace(1e-10,1,20000)])); sq=zz-q; Lq=AA/q+BB-3*q; Lp=AA/q**2+3; muq=q*(q-1)*(tt-q)
    J=ct(2*sq*Lq/muq,q,initial=0); y=np.exp(-J)*(1/w0+ct(np.exp(J)/muq,q,initial=0)); M=np.minimum(1/y,max(w0,2*Sh*(AA/zh+BB-3*zh)))
    I=ct((2*sq*Lq-M)/muq,q,initial=0); v=np.exp(I)*(w0-ct(2*sq*Lp*np.exp(-I),q,initial=0))
    ok+=v.min()>0
print('P2* corner',ok,len(out))

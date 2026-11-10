import numpy as np, pickle
bad=pickle.load(open('bad146.pkl','rb'))
def ok(P):
    tt,zz,AA,BB,Sh=P; zh=zz-Sh; mu_=zh*(zh-1)*(tt-zh); mup=(zh-1)*(tt-zh)+zh*(tt-zh)-zh*(zh-1); H=2*mu_+Sh*mup; w0=H/Sh
    Lh=AA/zh+BB-3*zh; th0=w0/(Sh*Lh)*0.999
    u=np.linspace(0,1,4001)[:-1]; q=zh+Sh*u; s=zz-q; L=AA/q+BB-3*q; Lp=AA/q**2+3; mu=q*(q-1)*(tt-q)
    for th1 in np.linspace(0.8,1.9,23):
        for u1 in (0.02,0.05,0.1,0.15,0.2,0.3,0.4,0.6):
            for th2 in np.linspace(0.2,th1,8):
                th=np.where(u<u1, th0+(th1-th0)*u/u1, th1+(th2-th1)*(u-u1)/(1-u1))
                thz=np.where(u<u1,(th1-th0)/(u1*Sh),(th2-th1)/((1-u1)*Sh))
                g=th*(2-th)*s*s*L*L+th*mu*L-(2-th)*mu*s*Lp-mu*thz*s*L
                if g.min()>0: return True
    return False
print(sum(ok(P) for P in bad),len(bad))

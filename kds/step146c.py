exec(open('step143.py').read().split("def two(P,th):")[0])
grid=[(a,b) for a in np.linspace(0.4,1.95,16) for b in np.linspace(0.1,1.9,19)]
def crit(P):
    tt,zz,AA,BB,Sh=P; zh=zz-Sh; mu_=zh*(zh-1)*(tt-zh); mup=(zh-1)*(tt-zh)+zh*(tt-zh)-zh*(zh-1); H=2*mu_+Sh*mup; w0=H/Sh
    Lh=AA/zh+BB-3*zh
    u=np.linspace(0,1,2001)[:-1]; q=zh+Sh*u; s=zz-q; L=AA/q+BB-3*q; Lp=AA/q**2+3; mu=q*(q-1)*(tt-q)
    for a,b in grid:
        if w0<a*Sh*Lh: continue
        th=a+(b-a)*u; thz=(b-a)/Sh
        if (th*(2-th)*s*s*L*L+th*mu*L-(2-th)*mu*s*Lp-mu*thz*s*L).min()>0: return (a,b)
    return None
bad=[Q for Q in out if crit(Q) is None]
print('corner fails',len(bad))
import pickle; pickle.dump(bad,open('bad146.pkl','wb'))

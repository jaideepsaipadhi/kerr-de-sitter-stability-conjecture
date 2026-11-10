exec(open('step143.py').read().split("def two(P,th):")[0])
from scipy.integrate import cumulative_trapezoid as ct, trapezoid
def crit(P,th=0.9):
    tt,zz,AA,BB,Sh=P; zh=zz-Sh; mu_=zh*(zh-1)*(tt-zh); mup=(zh-1)*(tt-zh)+zh*(tt-zh)-zh*(zh-1); H=2*mu_+Sh*mup; w0=H/Sh
    Lh=AA/zh+BB-3*zh; Phi0=H/Sh-2*Sh*Lh   # f(zhat)=1
    if Phi0>0: return True
    u=np.linspace(0,1,6001); q=zh+Sh*u; s=zz-q; L=AA/q+BB-3*q; Lp=AA/q**2+3; mu=q*(q-1)*(tt-q)
    D=s*s*L*L-2*mu*s*Lp
    k=np.minimum(w0,s*L+th*np.sqrt(np.maximum(D,0)))
    ok=(D>0)&(k< s*L+np.sqrt(np.maximum(D,0))+1e-300)
    # zone1 = initial run where ok and k nonincreasing
    j=len(u)-1
    bad=np.where(~ok)[0]
    if len(bad): j=bad.min()-1
    if j<2: return False
    kk=k[:j+1]
    mon=np.where(np.diff(kk)>1e-12*np.abs(kk[1:]))[0]
    if len(mon): j=mon.min()
    f=np.exp(ct(k[:j+1]/mu[:j+1],q[:j+1],initial=0))
    return Phi0+2*trapezoid(L[:j+1]*f,q[:j+1])>0
G=gen(1500,5)
for th in (0.0,0.5,0.9):
    print(th,'corner',sum(crit(Q,th) for Q in out),len(out),'generic',sum(crit(Q,th) for Q in G),len(G),flush=True)

exec(open('step143.py').read().split("def two(P,th):")[0])
ths=np.linspace(0.2,1.95,36)
def crit(P):
    tt,zz,AA,BB,Sh=P; zh=zz-Sh; mu_=zh*(zh-1)*(tt-zh); mup=(zh-1)*(tt-zh)+zh*(tt-zh)-zh*(zh-1); H=2*mu_+Sh*mup; w0=H/Sh
    Lh=AA/zh+BB-3*zh
    u=np.linspace(0,1,3001)[:-1]; q=zh+Sh*u; s=zz-q; L=AA/q+BB-3*q; Lp=AA/q**2+3; mu=q*(q-1)*(tt-q)
    for th in ths:
        if w0<th*Sh*Lh: continue
        if (th*(2-th)*s*s*L*L+th*mu*L-(2-th)*mu*s*Lp).min()>0: return th
    return None
G=gen(2000,5)
for nm,S in (('corner',out),('generic',G)):
    r=[crit(Q) for Q in S]; print(nm,len(S),'ok',sum(x is not None for x in r))

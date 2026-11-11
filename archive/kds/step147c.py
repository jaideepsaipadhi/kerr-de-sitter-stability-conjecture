exec(open('step143.py').read().split("def two(P,th):")[0])
ths=np.linspace(0.3,1.95,34)
def chk(P):
    tt,zz,AA,BB,Sh=P; zh=zz-Sh; mu_=zh*(zh-1)*(tt-zh); mup=(zh-1)*(tt-zh)+zh*(tt-zh)-zh*(zh-1); H=2*mu_+Sh*mup; w0=H/Sh
    Lh=AA/zh+BB-3*zh; ratio=w0/(Sh*Lh)
    x0=zz-1; q=1+x0*np.linspace(1e-6,1-1e-6,600); s=zz-q; L=AA/q+BB-3*q; Lp=AA/q**2+3; mu=q*(q-1)*(tt-q)
    for th in ths:
        if th<=ratio and (th*(2-th)*s*s*L*L+th*mu*L-(2-th)*mu*s*Lp).min()>0: return True
    return False
G=gen(2000,5)
print('corner',sum(chk(Q) for Q in out),len(out),'generic',sum(chk(Q) for Q in G),len(G))

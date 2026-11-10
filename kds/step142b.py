exec(open('step141b.py').read().split("for c in (0.5")[0])
def condL(P,th):
    tt,zz,AA,BB,Sh=P; zh=zz-Sh; mu_=zh*(zh-1)*(tt-zh); mup=(zh-1)*(tt-zh)+zh*(tt-zh)-zh*(zh-1); w0=(2*mu_+Sh*mup)/Sh
    if th*2*Sh*(AA/zh+BB-3*zh)>w0: return False
    sg=np.linspace(1e-6,1,1500); q=zz-Sh*sg; s=Sh*sg; L=AA/q+BB-3*q; Lp=AA/q**2+3; mu=q*(q-1)*(tt-q)
    return (-th*(1-th)*4*s**2*L**2-2*th*mu*(L+s*Lp)+2*mu*s*Lp).max()<=0
for th in (0.3,0.5,0.7,0.85,0.9):
    print(th,sum(condL(Q,th) for Q in out),len(out),flush=True)

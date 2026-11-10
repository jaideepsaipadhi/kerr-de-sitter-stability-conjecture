exec(open('step143.py').read().split("def two(P,th):")[0])
def crit(P):
    tt,zz,AA,BB,Sh=P; zh=zz-Sh; mu_=zh*(zh-1)*(tt-zh); mup=(zh-1)*(tt-zh)+zh*(tt-zh)-zh*(zh-1); H=2*mu_+Sh*mup
    Lh=AA/zh+BB-3*zh
    c1=H-Sh*Sh*Lh>0
    u=np.linspace(0,1,4001); q=zh+Sh*u; s=zz-q; L=AA/q+BB-3*q; Lp=AA/q**2+3; mu=q*(q-1)*(tt-q)
    C=s*s*L*L+mu*(L-s*Lp)
    return c1, C.min()>0
G=gen(2000,5)
for nm,S in (('corner',out),('generic',G)):
    r=np.array([crit(Q) for Q in S]); print(nm,len(S),'(i)',r[:,0].sum(),'(ii)',r[:,1].sum(),'both',(r[:,0]&r[:,1]).sum())

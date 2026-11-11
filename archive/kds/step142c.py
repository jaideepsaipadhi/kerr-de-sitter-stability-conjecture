exec(open('step141b.py').read().split("for c in (0.5")[0])
def condmax(P,c,th):
    tt,zz,AA,BB,Sh=P; zh=zz-Sh; mu_=zh*(zh-1)*(tt-zh); mup=(zh-1)*(tt-zh)+zh*(tt-zh)-zh*(zh-1); w0=(2*mu_+Sh*mup)/Sh
    sg=np.linspace(1e-6,1,3000); q=zz-Sh*sg; s=Sh*sg; L=AA/q+BB-3*q; Lp=AA/q**2+3; mu=q*(q-1)*(tt-q)
    k1=w0*sg**2*(1-c*sg*(1-sg)); k1z=-w0*(2*sg*(1-c*sg*(1-sg))+sg**2*(-c*(1-2*sg)))/Sh
    k2=th*2*s*L; k2z=-th*2*(L+s*Lp)
    if k2[-1]>w0: k2=np.minimum(k2,0)  # disallow
    act=k1>=k2; k=np.where(act,k1,k2); kz=np.where(act,k1z,k2z)
    return (mu*kz+k**2-2*s*L*k+2*mu*s*Lp).max()<=0
best={}
for c in (1.0,1.5):
  for th in (0.3,0.5,0.7):
    best[(c,th)]=[condmax(Q,c,th) for Q in out]
for k,v in sorted(best.items(),key=lambda kv:-sum(kv[1])): print(k,sum(v))
print('any',np.any(np.array(list(best.values())),axis=0).sum())

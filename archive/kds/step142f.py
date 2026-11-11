exec(open('step141b.py').read().split("for c in (0.5")[0])
def cm(P,c,th,a,generic=False):
    tt,zz,AA,BB,Sh=P; zh=zz-Sh; mu_=zh*(zh-1)*(tt-zh); mup=(zh-1)*(tt-zh)+zh*(tt-zh)-zh*(zh-1); w0=(2*mu_+Sh*mup)/Sh
    sg=np.linspace(1e-6,1,1200); q=zz-Sh*sg; s=Sh*sg; L=AA/q+BB-3*q; Lp=AA/q**2+3; mu=q*(q-1)*(tt-q)
    k1=a*w0*sg**2*(1-c*sg*(1-sg)); k1z=-a*w0*(2*sg*(1-c*sg*(1-sg))+sg**2*(-c*(1-2*sg)))/Sh
    D=s**2*L**2-2*mu*s*Lp; real=D>0; kp=np.where(real,s*L+np.sqrt(np.abs(D)),0)
    kpz=np.gradient(kp,q); k2=th*kp; k2z=th*kpz
    if k2[-1]>w0: return False
    act=k1>=k2; k=np.where(act,k1,k2); kz=np.where(act,k1z,k2z)
    return (mu*kz+k**2-2*s*L*k+2*mu*s*Lp)[1:-1].max()<=0
best={}
for c in (1.0,1.5):
  for th in (0.8,0.9,0.95):
    for a in (0.6,0.8,1.0):
      best[(c,th,a)]=[cm(Q,c,th,a) for Q in out]
for k,v in sorted(best.items(),key=lambda kv:-sum(kv[1]))[:5]: print(k,sum(v))
print('any',np.any(np.array(list(best.values())),axis=0).sum())

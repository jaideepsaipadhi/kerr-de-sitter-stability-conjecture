import numpy as np
rng=np.random.default_rng(1)
def sample(n):
    out=[]
    while len(out)<n:
        t=1+10**rng.uniform(-2,5); zf=rng.uniform(0,1); z0=1+zf*(t-1)
        B=(5.5*t-1.25)*10**rng.uniform(0,3); A=B*10**rng.uniform(-2,6)
        if 2*(z0-1)*(A+B-3)/(t-1)<=1: continue
        y=rng.uniform(0,1)**3 if rng.uniform()<.5 else 1-rng.uniform(0,1)**3
        out.append((t,z0,A,B,1+y*(z0-1)))
    return np.array(out)
P=sample(300000); t,z0,A,B,z=P.T
s=z0-z; x=z-1; x0=z0-1; u=x/x0
mu=z*x*(t-z); mup=-3*z*z+2*(1+t)*z-t
L=A/z+B-3*z; Lp=-A/z**2-3; Lpp=2*A/z**3
g=-Lp/L; gp=g*g-Lpp/L
kap=(A/z0+B-3*z0)/(z0*(z0-1)*(t-z0)); w=8*kap*s*s; psi=w/(1+w); psip=-16*kap*s/(1+w)**2
G=1.7*g; Gp=1.7*gp
q=0.45*(1-u)/(x*z); qp=0.45*(-(1/x0)/(x*z)-(1-u)*(2*z-1)/(x*z)**2)
h=G*psi+q; hp=Gp*psi+G*psip+qp
def rep(name,D,terms):
    n=sum(np.abs(T) for T in terms); m=D/n; i=np.argmax(m)
    print(name,'max normalised',m[i],'at t,zf,A/B,y=',t[i],(z0[i]-1)/(t[i]-1),A[i]/B[i],(z[i]-1)/(z0[i]-1))
T=[0.5*mup*h,0.5*mu*hp,0.5*mu*h*h,-s*L*h,-s*Lp]
rep('D exact',sum(T),T)
Tc=[0.5*mup*h,0.5*mu*Gp*psi,0.5*mu*h*h,-s*L*h,-s*Lp]
rep('D drop G psi\', q\'',sum(Tc),Tc)
Tc2=[0.5*mup*h,0.5*mu*(Gp*psi+qp),0.5*mu*h*h,-s*L*h,-s*Lp]
rep('D drop G psi\' only',sum(Tc2),Tc2)
for W0 in (1,3,10,30):
    msk=w>=W0
    Tb=[0.5*mup*h,0.5*mu*Gp*psi,0.5*mu*h*h,-s*L*h,-s*Lp]
    D=sum(Tb); n=sum(np.abs(T) for T in Tb); m=(D/n)[msk]; i=np.argmax(m)
    print('bulk w>=',W0,'drop Gpsi\',q\': max',m[i],' frac pts',msk.mean())

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
msk=w>=10
G2=1.7*g; qq=0.45/(x*z)
lhs=0.5*mu*(G2+qq)*(1/z+1/x)+0.85*mu*g*g+0.5*mu*(G2+qq)**2
rhs=0.5455*s*L*g+s*L*q
m=((lhs-rhs)/(lhs+rhs))[msk]; i=np.argmax(m)
idx=np.where(msk)[0][i]
print('(**) max normalised',m[i],'t,zf,A/B,y',t[idx],(z0[idx]-1)/(t[idx]-1),A[idx]/B[idx],(z[idx]-1)/(z0[idx]-1))
hup=G2+q; hlo=1.7*(10/11)*g+q
t1=np.where(mup>0,0.5*mup*hup,0.5*mup*hlo)
lhs2=t1+0.85*mu*g*g+0.5*mu*hup**2
rhs2=0.5455*s*L*g+s*L*q
m=((lhs2-rhs2)/(np.abs(t1)+0.85*mu*g*g+0.5*mu*hup**2+rhs2))[msk]; i=np.argmax(m); idx=np.where(msk)[0][i]
print('attempt2 max',m[i],'t,zf,A/B,y',t[idx],(z0[idx]-1)/(t[idx]-1),A[idx]/B[idx],(z[idx]-1)/(z0[idx]-1))
cpm=np.where(mup>0,1.7,17/11)
a_l=0.5*mup*cpm+2.295*mu*g; a_r=0.5455*s*L
b_l=0.5*mup+1.7*mu*g+0.5*mu*q; b_r=s*L
for nm,l_,r_ in (('(a)',a_l,a_r),('(b)',b_l,b_r)):
    m=((l_-r_)/(np.abs(l_)+r_))[msk]; i=np.argmax(m); idx=np.where(msk)[0][i]
    print(nm,'max',m[i],'t,zf,A/B,y',t[idx],(z0[idx]-1)/(t[idx]-1),A[idx]/B[idx],(z[idx]-1)/(z0[idx]-1))
for W0 in (1,3,10):
    mk=w>=W0
    hlo2=1.7*(W0/(1+W0))*g+q
    t1=np.where(mup>0,0.5*mup*hup,0.5*mup*hlo2)
    cW=1.7*W0/(1+W0)-1
    l3=t1+0.85*mu*g*g+0.5*mu*hup**2+0.5*mu*qp
    r3=cW*s*L*g+s*L*q
    m=((l3-r3)/(np.abs(t1)+0.85*mu*g*g+0.5*mu*hup**2+np.abs(0.5*mu*qp)+r3))[mk]; i=np.argmax(m); idx=np.where(mk)[0][i]
    print('B3 W0',W0,'max',m[i],'t,zf,A/B,y',t[idx],(z0[idx]-1)/(t[idx]-1),A[idx]/B[idx],(z[idx]-1)/(z0[idx]-1))
mk=w>=10; v=1-u
Iq_l=(t-z)*(0.225*(z+x)-0.12375)+0.765*mu*g; Iq_r=0.45*s*L
Ig_l=0.5*mup*cpm+2.295*mu*g; Ig_r=0.5455*s*L+(0.225*v+0.225*(t-z)/x0)/g
for nm,l_,r_ in (('Iq',Iq_l,Iq_r),('Ig',Ig_l,Ig_r)):
    m=((l_-r_)/(np.abs(l_)+np.abs(r_)))[mk]; i=np.argmax(m); idx=np.where(mk)[0][i]
    print(nm,'max',m[i],'t,zf,A/B,y',t[idx],(z0[idx]-1)/(t[idx]-1),A[idx]/B[idx],(z[idx]-1)/(z0[idx]-1))
N=A/z**2+3
l=0.85*(t-z)*(2*z-1)+2.295*z*x*(t-z)*N/L; r=0.5455*s*L+0.225*(t-z)*L/(x0*N)
m=((l-r)/(l+r))[mk]; i=np.argmax(m); idx=np.where(mk)[0][i]
print('Ig-crude max',m[i],'t,zf,A/B,y',t[idx],(z0[idx]-1)/(t[idx]-1),A[idx]/B[idx],(z[idx]-1)/(z0[idx]-1))
l=np.where(mup>0,0.85*mup,0.7727*mup)+2.295*z*x*(t-z)*N/L; r=0.5455*s*L+0.225*(t-z)*L/(x0*N)
m=((l-r)/(np.abs(l)+r))[mk]; i=np.argmax(m); idx=np.where(mk)[0][i]
print('Ig (exact mu\') max',m[i],'t,zf,A/B,y',t[idx],(z0[idx]-1)/(t[idx]-1),A[idx]/B[idx],(z[idx]-1)/(z0[idx]-1))
l=np.maximum(mup,0)*0.85+2.295*z*x*(t-z)*N/L
m=((l-r)/(np.abs(l)+r))[mk]; i=np.argmax(m); idx=np.where(mk)[0][i]
print('Ig (drop mu\'<0) max',m[i],'t,zf,A/B,y',t[idx],(z0[idx]-1)/(t[idx]-1),A[idx]/B[idx],(z[idx]-1)/(z0[idx]-1))
r2=0.5455*s*L
m=((l-r2)/(np.abs(l)+r2))[mk]; i=np.argmax(m); idx=np.where(mk)[0][i]
print('Ig (drop mu\'<0 and slack) max',m[i],'t,zf,A/B,y',t[idx],(z0[idx]-1)/(t[idx]-1),A[idx]/B[idx],(z[idx]-1)/(z0[idx]-1))
reg=mk&(z<=B/6)
lG=(t-z)*(3.995*z-3.145-0.225*z/x0)-0.7727*z*x; rG=0.5455*s*L
m=((lG-rG)/(np.abs(lG)+rG))[reg]; i=np.argmax(m); idx=np.where(reg)[0][i]
print('G1 max',m[i],'frac',reg.mean(),'t,zf,A/B,y',t[idx],(z0[idx]-1)/(t[idx]-1),A[idx]/B[idx],(z[idx]-1)/(z0[idx]-1))

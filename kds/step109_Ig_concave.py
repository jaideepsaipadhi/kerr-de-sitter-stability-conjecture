import numpy as np
rng=np.random.default_rng(5); worst1=worst2=-9; a1=a2=None; n1=n2=0
for _ in range(300000):
    t=1+10**rng.uniform(-2,5); z0=1+rng.uniform(0,1)*(t-1); z=1+rng.uniform(0,1)*(z0-1)
    B=5.5*t-1.25
    if z>=B/6: continue
    x=z-1; x0=z0-1; s=z0-z; mu=z*x*(t-z); mup=-3*z*z+2*(1+t)*z-t; mu0=z0*x0*(t-z0)
    P=np.where(mup>0,0.85,17/22)*mup
    b=0.5455*s; c=0.225*z*(t-z)/x0; d=B-6*z; a=2.295*mu*d/z
    K=P+2.295*mu/z-c
    # L_min: A>=0 -> L>=B-3z ; R3: 2x0 L(1)>t-1, 2x0^2(A/z0+3)>=t-1 ; bulk 8 L0 s^2>=10 mu0.  A = z(L-B+3z)
    # L(1)=A+B-3, L0=A/z0+B-3z0
    Amin=0.0
    Amin=max(Amin,(t-1)/(2*x0)-B+3+1e-12)
    Amin=max(Amin,z0*((t-1)/(2*x0*x0)-3))
    Amin=max(Amin,z0*(10*mu0/(8*s*s)-B+3*z0))
    Lmin=Amin/z+B-3*z
    F=lambda L: P+2.295*mu*(1/z-d/(z*L))-b*L-0.225*(t-z)/(x0*(1/z-d/(z*L)))
    Ls=d+np.sqrt(c*d/b)
    if Ls>=Lmin:
        v=(K-2*np.sqrt(b*c*d)-b*d)/(abs(K)+2*np.sqrt(b*c*d)+b*d); n1+=1
        if v>worst1: worst1=v; a1=(t,(z0-1)/(t-1),(z-1)/x0)
    v2=F(Lmin)/(abs(P)+2.295*mu/z+b*Lmin+c*Lmin/(Lmin-d)); n2+=1
    if v2>worst2: worst2=v2; a2=(t,(z0-1)/(t-1),(z-1)/x0)
print('I1 (interior, a/L dropped):',n1,worst1,a1); print('I2 (boundary L_min):',n2,worst2,a2)

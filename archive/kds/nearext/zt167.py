import numpy as np, time, sympy as sp, warnings; warnings.filterwarnings('ignore')
exec(open('bmin165.py').read().split("rng=np.random.default_rng(3)")[0].replace("NEGB1=[0]",""))
def geo(r0,r2,frac):
    r1=1.0; r3=-(r0+r1+r2); rs=[r0,r1,r2,r3]
    L,a,al,w1=kds(r0,r1,r2); Xi=1+al; a2=a*a
    Dp=lambda r: -(L/3)*sum(np.prod([r-x for k,x in enumerate(rs) if k!=j]) for j in range(4))
    P=[(x*x+a2)/abs(Dp(x)) for x in rs]; Dabs=[abs(Dp(x)) for x in rs]
    w2=a*(1/Dabs[1]-1/Dabs[0])/(P[1]-P[0]); w=w2+frac*(w1-w2)
    kap=[1j*Xi*(w*(x*x+a2)-a)/Dp(x) for x in rs]   # m=1
    t=(r1-r3)*(r2-r0)/((r1-r0)*(r2-r3)); X=[abs(kap[j]+kap[3]) for j in (0,1,2)]
    Ah=t*(X[0]**2-(X[1]+X[2])**2); C0=3/(L*(r1-r0)*(r2-r3)); c1=1+t+3*(1-al)*C0
    b2=C0*(1-2*Xi*Xi*a*w); b1=C0*(1-4*al*Xi-4*Xi*Xi*a*w); b0=C0*(-4-4*al)+c1
    e1=2*(kap[1]+kap[3]); et=2*(kap[2]+kap[3]); q1v=(t-1)*e1.imag/2; qtv=-t*(t-1)*et.imag/2
    z0=1+(t-1)*q1v/(q1v-qtv)
    return t,z0,Ah,b2,b1,b0
MV=[]; rng=np.random.default_rng(6); t0=time.time(); n=0; bad1=bad2=0; badg=0
while time.time()-t0<25:
    r0=1-10**rng.uniform(-3.5,0)*0.999; r2=1+10**rng.uniform(-4,3.5); frac=rng.uniform(0.0001,0.99999)
    try: tt,zz,Ah,b2,b1,b0=geo(r0,r2,frac)
    except Exception: continue
    x0=zz-1; zT=1+x0*np.linspace(1e-4,1-1e-4,300)
    A=AT(tt,zz,zT); m2=A/Ah; rec=np.concatenate([[True],A[1:]>np.maximum.accumulate(A)[:-1]])
    for k,Bl in ((1,m2*(b2+b1)+b0),(2,m2*b2+b0)):
        Nv=Nf(tt,zz,A,Bl,zT)
        R3=(2*x0*(A+Bl-3)>tt-1)&(2*x0**2*(A/zz+3)>=tt-1)&(A/zz+Bl-3*zz>0); g=(m2>=9)&R3&rec
        if k==1:
            bad1+=np.any((Nv>=0)&g); badg+=np.any(Nv>=0)
            if np.any((Nv>=0)&g): j=np.where((Nv>=0)&g)[0]; MV.append((m2[j].min(),m2[j].max(),(zT[j[0]]-1)/x0,b1))
        else: bad2+=np.any((Nv>=0)&g)
    n+=1
print(n,'viol Blow1 (m>=1)',bad1,' Blow2',bad2,' Blow1 ungated',badg)

MV=np.array(MV); print('m2 min over viol',np.percentile(MV[:,0],[0,50,100]),'b1>0 frac',(MV[:,3]>0).mean())

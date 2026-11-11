import flint
exec(open('smu_build.py').read().split("ex={'Q2'")[0])
FF=(f1,f2,Xi,lamv)
import numpy as np, warnings, os; warnings.filterwarnings('ignore')
os.chdir('../nearext'); exec(open('cos176.py').read().split("rng=np.random.default_rng(2)")[0])
def ex(r0,r1,r2,m,frac):
    r3=-(r0+r1+r2); rs=[r0,r1,r2,r3]
    L,a,al,w1=kds(r0,r1,r2); Xi=1+al; a2=a*a
    Dp=lambda r: -(L/3)*sum(np.prod([r-x for k,x in enumerate(rs) if k!=j]) for j in range(4))
    P=[(x*x+a2)/abs(Dp(x)) for x in rs]; Dabs=[abs(Dp(x)) for x in rs]
    w2=a*(1/Dabs[1]-1/Dabs[0])/(P[1]-P[0]); w=w2+frac*(w1-w2); om=m*w
    kap=[1j*Xi*(om*(x*x+a2)-a*m)/Dp(x) for x in rs]
    return (2*(kap[1]+kap[3])).imag,(2*(kap[2]+kap[3])).imag,a
for (r0v,r2v,fr) in [(0.5,3.0,0.3),(0.9,1.5,0.7),(0.2,20.0,0.5)]:
    u,et,av=ex(r0v,1.0,r2v,1,fr)
    pe=lambda p: sum(float(int(c))*r0v**int(e[0])*r2v**int(e[1])*(av*av)**int(e[2])*fr**int(e[3]) for e,c in p.to_dict().items())
    ev=lambda Fo: pe(Fo.n)/pe(Fo.d)
    f1v,f2v,Xv,lv=[ev(z) for z in FF]
    print(u/(Xv*av*f1v/lv), et/(Xv*av*f2v/lv))

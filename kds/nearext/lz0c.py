import numpy as np
def geom(rho,v):
    r0,r1,r2=rho,1.0,1/v; r3=-(r0+r1+r2); rs=[r0,r1,r2,r3]
    e2=sum(rs[i]*rs[j] for i in range(4) for j in range(i+1,4)); e4=r0*r1*r2*r3
    a_,b_,c_=e4/3,e2,3.0; L=(-b_-np.sqrt(b_*b_-4*a_*c_))/(2*a_)
    if L<=0: L=(-b_+np.sqrt(b_*b_-4*a_*c_))/(2*a_)
    a2=-L*e4/3; a=np.sqrt(a2); al=L*a2/3; Xi=1+al
    Dp=[-(L/3)*np.prod([x-y for k,y in enumerate(rs) if k!=j]) for j,x in enumerate(rs)]
    rho_=[x*x+a2 for x in rs]; P=[rho_[j]/abs(Dp[j]) for j in range(4)]; Om=[a/x for x in rho_]
    D=[abs(x) for x in Dp]
    w1=(P[1]*Om[1]+P[3]*Om[3])/(P[1]+P[3]); w2=a*(1/D[1]-1/D[0])/(P[1]-P[0])
    t=(r1-r3)*(r2-r0)/((r1-r0)*(r2-r3)); C0=3/(L*(r1-r0)*(r2-r3))
    return dict(rs=rs,L=L,a=a,al=al,Xi=Xi,P=P,Om=Om,w1=w1,w2=w2,t=t,C0=C0)
def coeffs(g,frac):
    P,Om,a,al,Xi,t,C0=g['P'],g['Om'],g['a'],g['al'],g['Xi'],g['t'],g['C0']
    w=g['w2']+frac*(g['w1']-g['w2'])
    u0=-P[0]*(w-Om[0]); u3=P[3]*(w-Om[3]); f1=P[1]*(w-Om[1])+u3; f2=-P[2]*(w-Om[2])+u3
    z0=t*(f1+f2)/(f1+t*f2)
    A2=Xi**2*t*((u0+u3)**2-(f1+f2)**2)
    Q2=A2/z0+C0*(1-2*Xi**2*a*w); Q1=C0*(1-4*al*Xi-4*Xi**2*a*w); Q0=-C0*(1+7*al)+1+t-3*z0
    return Q2,Q1,Q0,z0,A2,C0,t
def F(g,frac,m):
    Q2,Q1,Q0,z0,A2,C0,t=coeffs(g,frac)
    val=m*m*Q2+m*Q1+Q0+(4*C0 if m==1 else 0)
    return val/(m*m*abs(A2)/z0+m*m*C0+t)
for rho in [0.0,0.5,0.9,0.999,0.99999]:
  for v in [1e-5,1e-2,0.3,0.9,0.999]:
    try: g=geom(rho,v)
    except Exception as e: print('x',rho,v,e); continue
    vals=[min(F(g,fr,m) for fr in np.linspace(0,1,41)) for m in (1,2,3,10,1000)]
    print(rho,v,[round(x,3) for x in vals], 'Q2/norm', round(min(coeffs(g,fr)[0] for fr in np.linspace(0,1,41))/g['C0'],3))

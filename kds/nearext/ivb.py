import numpy as np, math, time
nx=np.nextafter
class I:
    __slots__=('lo','hi')
    def __init__(s,lo,hi=None):
        if hi is None: hi=lo
        s.lo=float(lo); s.hi=float(hi)
    def __add__(s,o): o=o if isinstance(o,I) else I(o); return I(nx(s.lo+o.lo,-np.inf),nx(s.hi+o.hi,np.inf))
    __radd__=__add__
    def __neg__(s): return I(-s.hi,-s.lo)
    def __sub__(s,o): o=o if isinstance(o,I) else I(o); return s+(-o)
    def __rsub__(s,o): return I(o)-s
    def __mul__(s,o):
        o=o if isinstance(o,I) else I(o); p=[s.lo*o.lo,s.lo*o.hi,s.hi*o.lo,s.hi*o.hi]
        if any(math.isnan(x) for x in p): return I(-np.inf,np.inf)
        return I(nx(min(p),-np.inf),nx(max(p),np.inf))
    __rmul__=__mul__
    def inv(s):
        if s.lo<=0<=s.hi: return I(-np.inf,np.inf)
        return I(nx(1/s.hi,-np.inf),nx(1/s.lo,np.inf))
    def __truediv__(s,o): o=o if isinstance(o,I) else I(o); return s*o.inv()
    def __rtruediv__(s,o): return I(o)*s.inv()
    def sqrt(s):
        if s.lo<0: return I(-np.inf,np.inf)
        return I(nx(math.sqrt(s.lo),-np.inf) if s.lo>0 else 0.0,nx(math.sqrt(s.hi),np.inf))
    def sq(s):
        if s.lo>=0: return s*s
        if s.hi<=0: return (-s)*(-s)
        m=max(-s.lo,s.hi); return I(0,nx(m*m,np.inf))
    def __abs__(s):
        if s.lo>=0: return s
        if s.hi<=0: return -s
        return I(0,max(-s.lo,s.hi))
def F(rho,v,fr,m):
    r0,r1=rho,I(1.0); r2=I(1.0)/v; r3=-(r0+r1+r2); rs=[r0,r1,r2,r3]
    e2=r0*r1+r0*r2+r0*r3+r1*r2+r1*r3+r2*r3; e4=r0*r1*r2*r3
    # Lambda positive root of (e4/3)L^2+e2 L+3=0, e4<0: L = (-e2 - sqrt(e2^2-4e4))/(2e4/3)  (positive)
    L=(-e2-(e2.sq()-(4/3*3)*e4).sqrt())/((2/3)*e4) if False else (e2+(e2.sq()-4*e4).sqrt())*(-1.5)/e4
    a2=-L*e4/3; a=a2.sqrt(); al=L*a2/3; Xi=1+al
    def prod(lst):
        p=I(1.0)
        for x in lst: p=p*x
        return p
    Dp=[-(L/3)*prod([rs[j]-rs[k] for k in range(4) if k!=j]) for j in range(4)]
    rr=[x.sq()+a2 for x in rs]; P=[rr[j]/abs(Dp[j]) for j in range(4)]; Om=[a/x for x in rr]; D=[abs(x) for x in Dp]
    w1=(P[1]*Om[1]+P[3]*Om[3])/(P[1]+P[3]); w2=a*(1/D[1]-1/D[0])/(P[1]-P[0])
    w=w2+fr*(w1-w2)
    t=(r1-r3)*(r2-r0)/((r1-r0)*(r2-r3)); C0=3/(L*(r1-r0)*(r2-r3))
    u0=-P[0]*(w-Om[0]); u3=P[3]*(w-Om[3]); f1=P[1]*(w-Om[1])+u3; f2=-P[2]*(w-Om[2])+u3
    z0=t*(f1+f2)/(f1+t*f2)
    A2=Xi.sq()*t*((u0+u3).sq()-(f1+f2).sq())
    Q2=A2/z0+C0*(1-2*Xi.sq()*a*w); Q1=C0*(1-4*al*Xi-4*Xi.sq()*a*w); Q0=-C0*(1+7*al)+1+t-3*z0
    val=m*m*Q2+m*Q1+Q0+(4*C0 if m==1 else 0)
    return val/C0, Q2/C0, (2*m*Q2+Q1)/C0
def bb(box,m,depth=0,cnt=[0]):
    cnt[0]+=1
    (a0,b0),(a1,b1),(a2,b2)=box
    val,q2,dq=F(I(a0,b0),I(a1,b1),I(a2,b2),m)
    ok = val.lo>0 and (m<3 or (q2.lo>0 and dq.lo>0))
    if ok: return True
    if depth>26: return False
    widths=[b0-a0,(b1-a1),(b2-a2)]; k=int(np.argmax(widths))
    mid=(box[k][0]+box[k][1])/2
    b1_=list(box); b2_=list(box); b1_[k]=(box[k][0],mid); b2_[k]=(mid,box[k][1])
    return bb(b1_,m,depth+1,cnt) and bb(b2_,m,depth+1,cnt)
if __name__=='__main__':
    t0=time.time()
    for m in (1,2,3):
        cnt=[0]; ok=bb([(0.05,0.95),(0.05,0.95),(0.0,1.0)],m,0,cnt); print(m,ok,cnt[0],round(time.time()-t0,1),flush=True)

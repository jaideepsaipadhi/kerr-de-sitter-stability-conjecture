import flint, time, pickle, numpy as np, sys
src=open('win_flint.py').read()
exec(src.split("print('R z0'")[0])
Xi=1-F(b*b,e4); al=F(-b*b,e4); lamv=F(-b,e4); K=F(one,(r1-r0)*(r2-r3))
Z=t/R
def P(x,k):
    r=F(one)
    for _ in range(k): r=r*x
    return r
c=lambda n,d=1: F(n*one,d*one)
NA=c(16,125)*Z*Z*(c(-4)*t*Z+c(5)*t+c(3)*Z*Z-c(4)*Z)
NB=c(16,3125)*P(Z,3)*(c(-75)*t*Z+c(100)*t+c(54)*Z*Z-c(75)*Z)
N0=c(2,3125)*Z*Z*(c(3375)*t*t*Z*Z-c(9000)*t*t*Z+c(5875)*t*t-c(3672)*t*P(Z,3)+c(11700)*t*Z*Z-c(9000)*t*Z+c(891)*P(Z,4)-c(3672)*P(Z,3)+c(3375)*Z*Z)
Am2=c(4)*Xi*Xi*F(b)*U[0]*U[3]*t/lamv          # lamv*A/m^2
b2=K*(1-2*Xi*Xi*F(b)*y); b1=K*(1-4*al*Xi-4*Xi*Xi*F(b)*y); b0=-K*(1+7*al)+lamv*(1+t)
q2=(-NA)*Am2+(-NB)*b2; q1=(-NB)*b1; q0=(-NB)*b0-lamv*N0
Hhalf=q2+q1*c(1,2)+q0*c(1,4)                   # m=2
H1=q2+q1+q0+(-NB)*(c(4)*K)                     # m=1 (l0=2)
gate=c(3)*Z-c(5)                               # 0.6 z0 > 1
out={k:(v.n,v.d) for k,v in dict(q2=q2,q1=q1,q0=q0,Hhalf=Hhalf,H1=H1,gate=gate,NB=NB,NA=NA).items()}
pickle.dump({k:(str(a),str(d)) for k,(a,d) in out.items()},open('zhat.pkl','wb'))
# numeric survey
rng=np.random.default_rng(0); res={k:[0,0] for k in out}
for _ in range(3000):
    x=1-10**rng.uniform(-3,0)*0.999; r2v=10**rng.uniform(0.05,3); tau=rng.uniform(0.001,0.999)
    r3v=-(x+1+r2v); e2v=x*1+x*r2v+x*r3v+r2v+r3v+r2v*r3v; e4v=x*r2v*r3v
    bv=(e2v+np.sqrt(e2v*e2v-4*e4v))/2
    vals={}
    for k,(a,d) in out.items():
        vals[k]=float(a(flint.fmpq(int(x*10**12),10**12),flint.fmpq(int(r2v*10**9),10**9),bv if False else flint.fmpq(0),flint.fmpq(0))) if False else None
    break
print('built', {k:(v[0].degrees(),len(v[0])) for k,v in out.items()})

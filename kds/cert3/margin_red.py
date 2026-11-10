import flint, time, pickle, numpy as np
t0=time.time()
src=open('win_flint.py').read()
exec(src.split("print('R z0'")[0])
Xi=1-F(b*b,e4); al=F(-b*b,e4); lamv=F(-b,e4); K=F(one,(r1-r0)*(r2-r3))
LQ2=2*Xi*Xi*F(b)*U[0]*U[3]*R/lamv + K*(1-2*Xi*Xi*F(b)*y)
LQ1=K*(1-4*al*Xi-4*Xi*Xi*F(b)*y)
LQ0=-K*(1+7*al)+lamv*(1+t-5*z0)
ex={'Q2':LQ2,'H2':LQ2+LQ1*F(1,2*one)+LQ0*F(1,4*one),'H1':LQ2+LQ1+LQ0+4*K}
def bcoeffs(p):
    d={}
    for exps,c in p.to_dict().items():
        k=exps[2]; e=(exps[0],exps[1],0,exps[3]); d.setdefault(k,{})[e]=c
    return {k:ctx.from_dict(v) for k,v in d.items()}
def reduce_b(p):
    cs=bcoeffs(p); deg=max(cs)
    for k in range(deg,1,-1):
        c=cs.pop(k,None)
        if c is None: continue
        # b^k = b^{k-2} * (e2 b - e4)
        cs[k-1]=cs.get(k-1,0*one)+c*e2; cs[k-2]=cs.get(k-2,0*one)-c*e4
    return cs.get(1,0*one), cs.get(0,0*one)
res={}
for k,e in ex.items():
    F1c,F0c=reduce_b(e.n*e.d)
    Y=-(2*F0c+F1c*e2); D=e2*e2-4*e4; P=F1c*F1c*D-Y*Y
    res[k]=(F1c,Y,P); print(k,'F1 deg',F1c.degrees(),'P deg',P.degrees(),'terms',len(P),round(time.time()-t0,1),flush=True)
pickle.dump({k:tuple(x.to_dict() for x in v) for k,v in res.items()},open('margin/win_res.pkl','wb'))
# numeric survey
def ev(p,a,bb,c):
    return float(p(flint.fmpz(0)+0,0,0,0)) if False else None
import random
for k,(F1c,Y,P) in res.items():
    s1=set(); sP=set(); sY=set()
    for _ in range(400):
        x=random.uniform(0.02,0.98); vv=random.uniform(0.02,0.98); tt=random.uniform(0,1)
        vals=[float(q(flint.fmpq(x).limit_denominator(10**6) if False else flint.fmpq(int(x*1e6),10**6), flint.fmpq(int(1e6/vv),10**6), 0, flint.fmpq(int(tt*1e6),10**6))) for q in (F1c,P,Y)]
        s1.add(np.sign(vals[0])); sP.add(np.sign(vals[1])); sY.add(np.sign(vals[2]))
    print(k,'F1',s1,'P',sP,'Y',sY, round(time.time()-t0,1))

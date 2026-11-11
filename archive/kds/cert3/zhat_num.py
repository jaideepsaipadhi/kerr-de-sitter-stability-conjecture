import flint, pickle, numpy as np, mpmath as mp
exec(open('zhat_build.py').read().split("out={k:")[0].replace("print","pass;0 and print"))
exprs=dict(q2=q2,q0=q0,Hhalf=Hhalf,H1=H1,gate=gate)
def mk(p):
    d=p.to_dict(); E=np.array([[int(x) for x in k] for k in d],float); C=[int(v) for v in d.values()]
    return E,C
num={k:(mk(v.n),mk(v.d)) for k,v in exprs.items()}
mp.mp.dps=50
def ev(EC,pt):
    E,C=EC; s=mp.mpf(0)
    for e,c in zip(E,C): s+=c*pt[0]**int(e[0])*pt[1]**int(e[1])*pt[2]**int(e[2])*pt[3]**int(e[3])
    return s
rng=np.random.default_rng(0); st={k:0 for k in exprs}; n=0; q0neg=0
import time; t0=time.time()
while time.time()-t0<45:
    x=mp.mpf(1-10**rng.uniform(-3,0)*0.999); r2v=mp.mpf(10**rng.uniform(0.05,3)); tau=mp.mpf(rng.uniform(0.001,0.999))
    r3v=-(x+1+r2v); e2v=x+x*r2v+x*r3v+r2v+r3v+r2v*r3v; e4v=x*r2v*r3v
    bv=(e2v+mp.sqrt(e2v*e2v-4*e4v))/2; pt=(x,r2v,bv,tau)
    v={k:ev(nn,pt)/ev(dd,pt) for k,(nn,dd) in num.items()}
    if v['gate']<=0: continue
    n+=1
    for k in ('q2','Hhalf','H1'): st[k]+=v[k]>0
    q0neg+=v['q0']<0
print('n',n,st,'q0<0',q0neg)

exec(open('step141b.py').read().split("for c in (0.5")[0])
import sys
def gen(n,seed):
    rng=np.random.default_rng(seed); o=[]
    while len(o)<n:
        tt=1+10**rng.uniform(-2,3); zz=1+(tt-1)*rng.uniform(1e-3,1-1e-3); x0=zz-1
        a0=max(0.0,zz*((tt-1)/(2*x0**2)-3)); AA=a0+10**rng.uniform(-3,4)*rng.integers(0,2)*(a0+1)
        bb=AA/(2*zz)+2*zz+10**rng.uniform(-4,1); BB=3*zz-AA/zz+bb
        if 2*x0*(AA+BB-3)<=tt-1: continue
        c=np.array([g(tt,zz,AA,BB) for g in cf],float); r=np.roots(c)
        r=[x.real for x in r if abs(x.imag)<1e-9*abs(x) and 0<x.real<x0]
        if len(r)!=1: continue
        o.append((tt,zz,AA,BB,r[0]))
    return o
def two(P,th):
    tt,zz,AA,BB,Sh=P; zh=zz-Sh; mu_=zh*(zh-1)*(tt-zh); mup=(zh-1)*(tt-zh)+zh*(tt-zh)-zh*(zh-1); w0=(2*mu_+Sh*mup)/Sh
    sg=np.linspace(1e-7,1,1500); q=zz-Sh*sg; s=Sh*sg; L=AA/q+BB-3*q; Lp=AA/q**2+3; mu=q*(q-1)*(tt-q)
    D=s*s*L*L-2*mu*s*Lp
    K1=s*L+th*np.sqrt(np.maximum(D,0))
    if K1[-1]>w0: return False
    if np.any(np.diff(K1[D>0])<0) and False: pass
    # zone1 valid where D>0 (strictly between roots, decreasing). try switch points i (tail on sg<=sg[i])
    cand=np.where(D>0)[0]
    # need zone1 valid on all of [i, end]: D>0 on indices >= i
    bad=np.where(D<=0)[0]; imin=(bad.max()+1) if len(bad) else 0
    for i in cand[cand>=imin][::max(1,len(cand)//60)]:
        ks=K1[i]; ss=s[:i+1]; mm=mu[:i+1]; LL=L[:i+1]; LLp=Lp[:i+1]
        rho=0.0
        for it in range(30):
            k=ks-rho*(s[i]-ss)
            if k.min()<=0: break
            need=((k*k-2*ss*LL*k+2*mm*ss*LLp)/mm).max()
            if need<=rho*(1+1e-9): return True
            rho=need*1.0000001

    return False
G=gen(1500,5)
for th in (0.5,0.8,0.95):
    print(th,'corner',sum(two(Q,th) for Q in out),len(out),'generic',sum(two(Q,th) for Q in G),len(G),flush=True)

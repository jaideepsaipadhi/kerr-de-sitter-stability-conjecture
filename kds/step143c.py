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
    if s[-1]*L[-1]>w0: return False
    # zone1 valid where D>0 (strictly between roots, decreasing). try switch points i (tail on sg<=sg[i])
    cand=np.where(D>0)[0]
    # need zone1 valid on all of [i, end]: D>0 on indices >= i
    bad=np.where(D<=0)[0]; imin=(bad.max()+1) if len(bad) else 0
    for i in cand[cand>=imin][::max(1,len(cand)//60)]:
        c=L[i]; ss=s[:i+1]
        k=c*ss; R=-k*k+2*ss*L[:i+1]*k-2*mu[:i+1]*ss*Lp[:i+1]
        if (mu[:i+1]*(-c)-R)[1:].max()<=0: return True
    return False
G=gen(1500,5)
print('corner',sum(two(Q,0) for Q in out),len(out),'generic',sum(two(Q,0) for Q in G),len(G),flush=True)

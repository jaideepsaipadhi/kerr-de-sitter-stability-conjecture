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
    sg=np.linspace(1e-7,1,3000); q=zz-Sh*sg; s=Sh*sg; L=AA/q+BB-3*q; Lp=AA/q**2+3; mu=q*(q-1)*(tt-q)
    L0=AA/zz+BB-3*zz
    D=s*s*L*L-2*mu*s*Lp
    k1=np.where(D>0,s*L,-np.inf); k2=th*L0*s
    k=np.maximum(k1,k2); kz=np.where(k1>=k2,-(L+s*Lp),-th*L0)
    if k[-1]>w0: return False
    R=-k*k+2*s*L*k-2*mu*s*Lp
    return (mu*kz-R)[1:].max()<=0
G=gen(1500,5)
for th in (0.5,1.0,1.5,2.0):
    print(th,'corner',sum(two(Q,th) for Q in out),len(out),'generic',sum(two(Q,th) for Q in G),len(G),flush=True)

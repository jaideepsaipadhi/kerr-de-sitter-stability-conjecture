import numpy as np, pickle, sys, flint, time, sympy as sp
exec(open('win_bern.py').read().split("which=")[0])
x,v,ta,mu,nu=sp.symbols('x v tau mu nu')
def tosym(A):
    return sum(int(A[i,j,k])*x**i*v**j*ta**k for i,j,k in zip(*np.nonzero(A)))
def strip_pow(p,var):
    P=sp.Poly(p,var); k=min(m[0] for m in P.monoms()); return sp.expand(sp.cancel(p/var**k)),k
def cert(expr,V,label,tmax=25):
    P=sp.Poly(sp.expand(expr),*V)
    d=P.as_dict(); n=[P.degree(z) for z in V]
    A=np.zeros([q+1 for q in n],dtype=object)
    for e,c in d.items(): A[e]=int(c)
    B=bern(A); ok,cl,info=bb(B,tmax=tmax,maxdep=30); print(label,ok,cl,info,round(time.time()-t0,1),flush=True); return ok
which=sys.argv[1]; nm=sys.argv[2]
F1,Y,Pd=res[which]; d=F1 if nm=='F1' else Pd
c,fl=ctx.from_dict(d).factor()
for f,m in fl:
    if sum(f.degrees())<=8: continue
    p=tosym(to_arr(f))
    # chart A: tau = 1 - nu, v = mu*nu ; chart B: tau = 1 - mu*nu, v = nu
    pa,ka=strip_pow(sp.expand(p.subs({ta:1-nu,v:mu*nu})),nu)
    pb,kb=strip_pow(sp.expand(p.subs({ta:1-mu*nu,v:nu})),nu)
    print(which,nm,f.degrees(),'strip',ka,kb,flush=True)
    cert(pa,(x,mu,nu),' chartA'); cert(pb,(x,mu,nu),' chartB')

import numpy as np
rng=np.random.default_rng(41); st={'n':0,'Nx0neg':0,'concave':0,'convex':0,'neither':0,'singleroot':0}
for _ in range(200000):
    t=1+10**rng.uniform(-2,5); zf=rng.uniform(0,1); z0=1+zf*(t-1); x0=z0-1
    B=(5.5*t-1.25)*10**rng.uniform(0,3); A=B*10**rng.uniform(-3,6)
    if 2*x0*(A+B-3)/(t-1)<=1 or 2*x0**2*(A/z0+3)<t-1: continue
    mu0=z0*x0*(t-z0); mup0=-3*z0**2+2*(1+t)*z0-t
    c=[6*mu0**2, -6*mu0*mup0,
       4*(z0-1)*(z0-t)*(A+B*z0)+6*t**2*z0**2-6*t**2*z0+t**2-8*t*z0**3+12*t*z0**2-6*t*z0+3*z0**4-8*z0**3+6*z0**2,
       2*(A*(t+1-2*z0)+B*mup0+3*z0*(-3*t*z0+2*t+4*z0**2-3*z0)), 0, 2*(B+4*t-15*z0+4), 15]
    P=np.polynomial.Polynomial(c); S=np.linspace(0,x0,3001)[1:-1]
    st['n']+=1
    if P(x0)<0: st['Nx0neg']+=1
    d2=P.deriv(2)(S)
    if np.all(d2<0): st['concave']+=1
    elif np.all(d2>0): st['convex']+=1
    else: st['neither']+=1
    v=P(S); st['singleroot']+= (np.sum(np.diff(np.sign(v))!=0)==1)
print(st)
rng=np.random.default_rng(42); bad=[]
for _ in range(100000):
    t=1+10**rng.uniform(-2,5); zf=rng.uniform(0,1); z0=1+zf*(t-1); x0=z0-1
    B=(5.5*t-1.25)*10**rng.uniform(0,3); A=B*10**rng.uniform(-3,6)
    if 2*x0*(A+B-3)/(t-1)<=1 or 2*x0**2*(A/z0+3)<t-1: continue
    mu0=z0*x0*(t-z0); mup0=-3*z0**2+2*(1+t)*z0-t
    c=[6*mu0**2, -6*mu0*mup0,
       4*(z0-1)*(z0-t)*(A+B*z0)+6*t**2*z0**2-6*t**2*z0+t**2-8*t*z0**3+12*t*z0**2-6*t*z0+3*z0**4-8*z0**3+6*z0**2,
       2*(A*(t+1-2*z0)+B*mup0+3*z0*(-3*t*z0+2*t+4*z0**2-3*z0)), 0, 2*(B+4*t-15*z0+4), 15]
    # scale S = x0*u
    cu=[c[k]*x0**k for k in range(7)]; m=max(abs(v) for v in cu); cu=[v/m for v in cu]
    r=np.polynomial.Polynomial(cu).roots(); rr=[x.real for x in r if abs(x.imag)<1e-9 and 0<x.real<1]
    if len(rr)!=1: bad.append((t,zf,A/B,B/(5.5*t),len(rr)))
print('bad',len(bad),bad[:8])
from numpy.polynomial import polynomial as Pl
rng=np.random.default_rng(43); counts={}
for _ in range(100000):
    t=1+10**rng.uniform(-2,5); zf=rng.uniform(0,1); z0=1+zf*(t-1); x0=z0-1
    B=(5.5*t-1.25)*10**rng.uniform(0,3); A=B*10**rng.uniform(-3,6)
    if 2*x0*(A+B-3)/(t-1)<=1 or 2*x0**2*(A/z0+3)<t-1: continue
    mu0=z0*x0*(t-z0); mup0=-3*z0**2+2*(1+t)*z0-t
    c=[6*mu0**2, -6*mu0*mup0,
       4*(z0-1)*(z0-t)*(A+B*z0)+6*t**2*z0**2-6*t**2*z0+t**2-8*t*z0**3+12*t*z0**2-6*t*z0+3*z0**4-8*z0**3+6*z0**2,
       2*(A*(t+1-2*z0)+B*mup0+3*z0*(-3*t*z0+2*t+4*z0**2-3*z0)), 0, 2*(B+4*t-15*z0+4), 15]
    cu=[c[k]*x0**k for k in range(7)]
    # u = v/(1+v):  sum cu_k v^k (1+v)^(6-k)
    tot=np.zeros(7)
    for k in range(7):
        tot[:]+=np.pad(cu[k]*Pl.polypow([1,1],6-k)*np.eye(1,7-k,0).ravel()[0] if False else 0,(0,0))
    tot=np.zeros(7)
    for k in range(7):
        p=Pl.polymul(np.eye(1,k+1,k).ravel(),Pl.polypow([1,1],6-k))*cu[k]
        tot[:len(p)]+=p
    sg=np.sign(tot[np.abs(tot)>1e-12*np.abs(tot).max()])
    ch=int(np.sum(sg[1:]!=sg[:-1])); counts[ch]=counts.get(ch,0)+1
print('sign changes of Mobius-transformed coefficients:',counts)
rng=np.random.default_rng(44); sig=np.zeros((7,3),int)
for _ in range(100000):
    t=1+10**rng.uniform(-2,5); zf=rng.uniform(0,1); z0=1+zf*(t-1); x0=z0-1
    B=(5.5*t-1.25)*10**rng.uniform(0,3); A=B*10**rng.uniform(-3,6)
    if 2*x0*(A+B-3)/(t-1)<=1 or 2*x0**2*(A/z0+3)<t-1: continue
    mu0=z0*x0*(t-z0); mup0=-3*z0**2+2*(1+t)*z0-t
    c=[6*mu0**2, -6*mu0*mup0,
       4*(z0-1)*(z0-t)*(A+B*z0)+6*t**2*z0**2-6*t**2*z0+t**2-8*t*z0**3+12*t*z0**2-6*t*z0+3*z0**4-8*z0**3+6*z0**2,
       2*(A*(t+1-2*z0)+B*mup0+3*z0*(-3*t*z0+2*t+4*z0**2-3*z0)), 0, 2*(B+4*t-15*z0+4), 15]
    cu=[c[k]*x0**k for k in range(7)]; tot=np.zeros(7)
    for k in range(7):
        p=Pl.polymul(np.eye(1,k+1,k).ravel(),Pl.polypow([1,1],6-k))*cu[k]; tot[:len(p)]+=p
    for j in range(7): sig[j,int(np.sign(tot[j]))+1]+=1
print('rows j: counts of (neg, zero, pos):'); print(sig)

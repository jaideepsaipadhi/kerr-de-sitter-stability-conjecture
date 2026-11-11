import flint, time, pickle, sys, sympy as sp
src=open('win_flint.py').read()
exec(src.split("print('R z0'")[0])
zs,ts,z0s,As,Bs=sp.symbols('z t z0 A B')
s_=z0s-zs; mus=zs*(zs-1)*(ts-zs); Ls=As/zs+Bs-3*zs; Fs=mus/s_**2
nu,dn=sp.fraction(sp.factor(sp.together(sp.diff(mus*sp.diff(Fs,zs),zs)-2*s_*sp.diff(Ls*Fs,zs))))
Nexp=sp.expand(nu)
parts={'NA':sp.Poly(Nexp.coeff(As),ts,z0s,zs),'NB':sp.Poly(Nexp.coeff(Bs),ts,z0s,zs),'N0':sp.Poly(Nexp.subs({As:0,Bs:0}),ts,z0s,zs)}
Xi=1-F(b*b,e4); al=F(-b*b,e4); lamv=F(-b,e4); K=F(one,(r1-r0)*(r2-r3))
Z=t/R
c=lambda n,d=1: F(n*one,d*one)
zst=1+c(31,50)*(Z-1)*(Z-1)/Z
def P(x,k):
    r=F(one)
    for _ in range(k): r=r*x
    return r
cache={}
def pw(x,k,key):
    if (key,k) not in cache: cache[(key,k)]=P(x,k)
    return cache[(key,k)]
def conv(poly):
    acc=None
    for (i,j,k),cf in zip(poly.monoms(),poly.coeffs()):
        term=c(int(cf))*pw(t,i,'t')*pw(Z,j,'Z')*pw(zst,k,'z')
        acc=term if acc is None else acc+term
    return acc
t0=time.time()
NA=conv(parts['NA']); NB=conv(parts['NB']); N0=conv(parts['N0'])
print('conv',round(time.time()-t0,1),flush=True)
Am2=c(4)*Xi*Xi*F(b)*U[0]*U[3]*t/lamv
b2=K*(1-2*Xi*Xi*F(b)*y); b1=K*(1-4*al*Xi-4*Xi*Xi*F(b)*y); b0=-K*(1+7*al)+lamv*(1+t)
q2=(-NA)*Am2+(-NB)*b2; q1=(-NB)*b1; q0=(-NB)*b0-lamv*N0
Hhalf=q2+q1*c(1,2)+q0*c(1,4); H1=q2+q1+q0+(-NB)*(c(4)*K)
print('built',round(time.time()-t0,1),{k:v.n.degrees() for k,v in dict(q2=q2,q0=q0,Hhalf=Hhalf,H1=H1).items()},flush=True)

exec(open("sal_build.py").read().split("\nex={")[0])
T=t-z0; aa=z0-1
Gam=F(4*one)*U[0]*U[3]*T*z0*z0/(t*(f1+f2)*(f1+f2))
den=F(4*one)*Xi*Xi*F(b)*f2*(f1+t*f2)
Hm={1:LQ2p+LQ1+LQ0p+4*K, 2:LQ2p+LQ1*F(1,2*one)+LQ0p*F(1,4*one)}
Q={'gam':Gam/T,'alp':aa/T}
for m,H in Hm.items(): Q[f'rho{m}']=F(21,10*one)*(-H)*lamv/den
g=ctx.gens(); one_=one
for k,q in Q.items():
    n1=q.n.compose(one_,g[1],g[2],g[3],ctx=ctx); d1=q.d.compose(one_,g[1],g[2],g[3],ctx=ctx)
    print(k,'n deg',q.n.degrees(),'at r0=1: n zero?',n1==0,'d zero?',d1==0,flush=True)

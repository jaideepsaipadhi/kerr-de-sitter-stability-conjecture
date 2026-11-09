import sympy as sp, pickle, sys, time
T,ze,y,p,A,B=sp.symbols('T ze y p A B',positive=True); V=(T,ze,y,p)
gen=pickle.load(open('gen.pkl','rb'))
t=1/T; floor=sp.Rational(11,2)*t-sp.Rational(5,4); x0=ze*(t-1); z0=1+x0; s=(1-y)*x0; mu0=z0*x0*(t-z0)
SURF={
 'bulkB': lambda: (z0*(10*mu0/(8*s**2)-floor/p+3*z0), floor/p),
 'ppB':   lambda: ((t-1)/(2*x0)-floor/p+3, floor/p),
 'c2B':   lambda: (z0*((t-1)/(2*x0**2)-3), floor/p),
 'bulkw': lambda: (None,None),
}
def compose(P,an,ad,bn,bd):
    an,ad,bn,bd=[sp.Poly(sp.expand(e),*V) for e in (an,ad,bn,bd)]
    dA=P.degree(A); dB=P.degree(B)
    pa={i:an**i*ad**(dA-i) for i in range(dA+1)}; pb={j:bn**j*bd**(dB-j) for j in range(dB+1)}
    # group P by (i,j)
    from collections import defaultdict
    grp=defaultdict(lambda: sp.Poly(0,*V))
    for mon,co in zip(P.monoms(),P.coeffs()):
        i,j,k1,k2,k3=mon; grp[(i,j)]=grp[(i,j)]+sp.Poly(co*T**k1*ze**k2*y**k3,*V)
    out=sp.Poly(0,*V)
    for (i,j),c in grp.items(): out=out+c*pa[i]*pb[j]
    return out,dA,dB
def _old(P,an,ad,bn,bd):
    dA=P.degree(A); dB=P.degree(B); out=sp.Poly(0,*V)
    for (i,j,a,b,c),co in zip(P.monoms(),P.coeffs()):
        pass
    acc=0
    for mon,co in zip(P.monoms(),P.coeffs()):
        i,j,k1,k2,k3=mon
        acc+= co*T**k1*ze**k2*y**k3*an**i*ad**(dA-i)*bn**j*bd**(dB-j)
    return sp.Poly(sp.expand(acc),*V),dA,dB
name=sys.argv[1]; t0=time.time()
Aex,Bex=SURF[name]()
an,ad=sp.fraction(sp.together(Aex)); bn,bd=sp.fraction(sp.together(Bex))
res={}
for k,(P,dn) in gen.items():
    Q,dA,dB=compose(P,an,ad,bn,bd)
    den=dn.subs({A:Aex,B:Bex})*ad**dA*bd**dB
    res[k]=(Q,den); print(name,k,[Q.degree(v) for v in V],round(time.time()-t0,1),flush=True)
pickle.dump(res,open(f'surf_{name}.pkl','wb'))

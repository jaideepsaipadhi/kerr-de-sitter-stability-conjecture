import sympy as sp, pickle, numpy as np
r0,r2,lam=sp.symbols('r0 r2 lam',positive=True)
raw=pickle.load(open('w2_raw.pkl','rb')); red,e2,e4=pickle.load(open('w2_red.pkl','rb'))
f=lambda e,v=(r0,r2): sp.lambdify(v,e,'numpy')
R0,R2=np.meshgrid(np.linspace(0.02,0.98,50),1/np.linspace(0.02,0.98,50))
E2,E4=f(e2)(R0,R2),f(e4)(R0,R2); s=np.sqrt(E2**2-4*E4); L=(-E2-s)/(2*E4)
for k,(n,d) in raw.items():
    kdeg=sp.Poly(n,lam).degree()-1
    nv=f(n,(r0,r2,lam))(R0,R2,L); dv=f(d,(r0,r2,lam))(R0,R2,L)
    F1,F0=red[k]; pv=f(F1)(R0,R2)*L+f(F0)(R0,R2)
    print(k,'deg_lam',kdeg+1,'expr>0 all',bool((nv/dv>0).all()),' sign(pr)/sign(n) consistent with e4^k:',bool((np.sign(pv)==np.sign(nv)*np.sign(E4)**kdeg).all()))

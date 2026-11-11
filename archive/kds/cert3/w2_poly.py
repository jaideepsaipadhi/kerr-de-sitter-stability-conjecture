import sympy as sp, pickle, time, numpy as np
t0=time.time()
r0,r2,lam,v=sp.symbols('r0 r2 lam v',positive=True)
red,e2,e4=pickle.load(open('w2_red.pkl','rb'))
D=sp.expand(e2**2-4*e4); out={}
for k,(F1,F0) in red.items():
    X=sp.expand(-F1*e2+2*e4*F0)
    P=sp.expand(F1**2*D-X**2)
    for nm,e in (('F1',F1),('P',P)):
        ev=sp.expand(e.subs(r2,1/v)); dg=sp.Poly(sp.expand(e),r2).degree()
        ev=sp.expand(sp.cancel(ev*v**dg))
        fac=sp.factor_list(ev)
        out[(k,nm)]=fac
        print(k,nm,'deg',sp.Poly(ev,r0,v).degree(r0),sp.Poly(ev,r0,v).degree(v),'factors',[ (sp.Poly(f,r0,v).total_degree(),m) for f,m in fac[1]], 'const',fac[0], round(time.time()-t0,1),flush=True)
pickle.dump(out,open('w2_poly.pkl','wb'))

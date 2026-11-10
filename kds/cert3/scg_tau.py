import flint, time, pickle
exec(open("scg_build.py").read().split("\nres={}")[0])
def tau_coeffs(p):
    d={}
    for e,c in p.to_dict().items():
        k=int(e[3]); d.setdefault(k,{})[(int(e[0]),int(e[1]),int(e[2]),0)]=c
    return {k:ctx.from_dict(v) for k,v in d.items()}
out={}
for k,e in ex.items():
    cd,fld=e.d.factor(); assert all(m%2==0 for f,m in fld),'d not square'
    N=e.n*(1 if cd>0 else -1); print(k,'d const',cd); cs=tau_coeffs(N); print(k,'tau-degree',max(cs),flush=True)
    c0=cs.get(0,0*one); c1=cs.get(1,0*one); c2=cs.get(2,0*one)
    conds={'q0':c0,'q1':c0+c1+c2,'c2neg':-c2,'disc':4*c0*c2-c1*c1,'v1':c1,'v2':-(c1+2*c2),'mm':c0-c1-c2,'pm':c0+c1-c2,'mp':c0-c1+c2}
    for nm,p in conds.items():
        F1c,F0c=reduce_b(p); Y=-(2*F0c+F1c*e2); D=e2*e2-4*e4; P=F1c*F1c*D-Y*Y
        out[(k,nm)]=(F1c.to_dict(),Y.to_dict(),P.to_dict())
        print(' ',nm,F1c.degrees(),P.degrees(),round(time.time()-t0,1),flush=True)
pickle.dump(out,open('margin/scg_tau.pkl','wb'))

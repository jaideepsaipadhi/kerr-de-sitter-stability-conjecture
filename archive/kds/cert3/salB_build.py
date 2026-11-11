import flint, time, pickle
exec(open('smu_build.py').read().split("ex={'Q2'")[0])
two=F(2*one); uh=two*Xi*f1/lamv; eh=two*Xi*f2/lamv; Fh=-(uh+t*eh)
G=two*Fh+(t-1)*(Fh-F(3*one)*uh)
P2=F(3*one)*t*(t-1)*F(b)*eh*eh*Fh*lamv
LQ2p=K*(1-2*Xi*Xi*F(b)*y)
c=F(21,10*one)*G
ex={'Q2':P2+c*LQ2p,'H2':P2+c*(LQ2p+LQ1*F(1,2*one)+LQ0p*F(1,4*one)),'H1':P2+c*(LQ2p+LQ1+LQ0p+4*K)}
exec("def bcoeffs"+open('margin_red.py').read().split("def bcoeffs")[1].split("res={}")[0])
res={}
for k,e in ex.items():
    cd,fld=e.d.factor(); cn,fln=e.n.factor(); N=one*(int(cn)*(1 if cd>0 else -1))
    for f,m in list(fln)+list(fld):
        if m%2==1: N=N*f
        else: print(k,'drop square',f.degrees(),m)
    F1c,F0c=reduce_b(N); Y=-(2*F0c+F1c*e2); D=e2*e2-4*e4; P=F1c*F1c*D-Y*Y
    res[k]=(F1c,Y,P); print(k,F1c.degrees(),P.degrees(),flush=True)
pickle.dump({k:tuple(x.to_dict() for x in v) for k,v in res.items()},open('margin/sal_res.pkl','wb'))

import flint, time, pickle
exec(open('smu_build.py').read().split("ex={'Q2'")[0])
Ql=F(4*one)*Xi*Xi*F(b)*(f1+t*f2)*(f1+t*f2)/lamv; LQ2p=K*(1-2*Xi*Xi*F(b)*y)
c=F(6*one)*t
ex={'Q2':Ql+c*LQ2p,'H2':Ql+c*(LQ2p+LQ1*F(1,2*one)+LQ0p*F(1,4*one)),'H1':Ql+c*(LQ2p+LQ1+LQ0p+4*K)}
exec("def bcoeffs"+open('margin_red.py').read().split("def bcoeffs")[1].split("res={}")[0])
res={}
for k,e in ex.items():
    F1c,F0c=reduce_b(e.n*e.d); Y=-(2*F0c+F1c*e2); D=e2*e2-4*e4; P=F1c*F1c*D-Y*Y
    res[k]=(F1c,Y,P); print(k,F1c.degrees(),P.degrees(),flush=True)
pickle.dump({k:tuple(x.to_dict() for x in v) for k,v in res.items()},open('margin/scg6_res.pkl','wb'))

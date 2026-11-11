import flint, time, pickle
exec(open('smu_build.py').read().split("ex={'Q2'")[0])
Q2l=t*t*Xi*Xi*F(b)*(f1+f2)*(f1+f2)/(lamv*z0*z0)   # = lamv*q1^2/m^2
T=t-z0; aa=z0-1; c=F(105,100*one)
def W(coef,Tp):  # coef*lamv*Q2*Tp - c*t(t-1)*(-lamv L)  >=0 ; per m-power
    base=F(coef*one)*Q2l*Tp
    return {'Q2':base+c*t*(t-1)*LQ2p,'H2':base+c*t*(t-1)*(LQ2p+LQ1*F(1,2*one)+LQ0p*F(1,4*one)),'H1':base+c*t*(t-1)*(LQ2p+LQ1+LQ0p+4*K)}
ex={}
for k,v in W(2,T).items(): ex['C'+k]=v
for k,v in W(6,T*T*(1/(F(2*one)+F(4*one)*aa+T))).items(): ex['A'+k]=v
exec("def bcoeffs"+open('margin_red.py').read().split("def bcoeffs")[1].split("res={}")[0])
res={}
for k,e in ex.items():
    F1c,F0c=reduce_b(e.n*e.d); Y=-(2*F0c+F1c*e2); D=e2*e2-4*e4; P=F1c*F1c*D-Y*Y
    res[k]=(F1c,Y,P); print(k,F1c.degrees(),P.degrees(),round(time.time()-t0,1),flush=True)
pickle.dump({k:tuple(x.to_dict() for x in v) for k,v in res.items()},open('margin/sca_res.pkl','wb'))
pickle.dump({'Q2l':Q2l,'ex':ex} if False else None,open('/dev/null','wb'))

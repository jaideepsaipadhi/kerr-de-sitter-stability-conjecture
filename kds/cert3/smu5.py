import flint, time, pickle
t0=time.time()
src=open('win5.py').read(); exec(src.split("print('R z0'")[0])
exec(open('margin_red5.py').read().split("ex={")[0].split("exec(src.split")[1].split("\n",1)[1])
LQ2p=F(4*one)*Xi*Xi*F(b)*U[0]*U[3]/lamv + K*(1-2*Xi*Xi*F(b)*y)
LQ0p=-K*(1+7*al)+lamv*(1-2*t)
ex={'Q2':LQ2p,'H2':LQ2p+LQ1*F(1,2*one)+LQ0p*F(1,4*one),'H1':LQ2p+LQ1+LQ0p+4*K,'G':F(2*one)*(z0-1)-(t-1)}
exec("def bcoeffs"+open('margin_red5.py').read().split("def bcoeffs")[1].split("res={}")[0])
res={}
for k,e in ex.items():
    F1c,F0c=reduce_b(e.n*e.d)
    Y=-(2*F0c+F1c*e2); D=e2*e2-4*e4; P=F1c*F1c*D-Y*Y
    res[k]=(F1c,Y,P); print(k,'F1',F1c.degrees(),'P',P.degrees(),len(P),round(time.time()-t0,1),flush=True)
pickle.dump({k:tuple(x.to_dict() for x in v) for k,v in res.items()},open('margin/smu_res.pkl','wb'))

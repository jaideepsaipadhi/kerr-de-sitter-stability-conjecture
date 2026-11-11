import pickle, mpmath as mp
mp.mp.dps=80
T=pickle.load(open('margin/sal_res.pkl','rb'))
def ev(d,r0,r2,ta): return mp.fsum(int(c)*r0**int(e[0])*r2**int(e[1])*ta**int(e[3]) for e,c in d.items())
for (x,v,ta) in [('0.9858','0.8643','0.9763'),('0.3464','0.0001','0.9998')]:
    r0=mp.mpf(x); r2=1/mp.mpf(v); ta=mp.mpf(ta); r3=-(r0+1+r2); rs=[r0,1,r2,r3]
    E2=mp.fsum(rs[i]*rs[j] for i in range(4) for j in range(i+1,4)); E4=r0*r2*r3; D=E2*E2-4*E4
    for k in ('Q2','H2'):
        F1,Y,P=[ev(T[k][i],r0,r2,ta) for i in range(3)]
        print(x,v,k,mp.nstr(F1,4),mp.nstr(Y,4),mp.nstr(P,4),'val',mp.nstr((F1*mp.sqrt(D)-Y)/2,4))

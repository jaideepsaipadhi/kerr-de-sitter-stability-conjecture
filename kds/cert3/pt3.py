import pickle, mpmath as mp
mp.mp.dps=60
R=pickle.load(open('margin/sbe2_res.pkl','rb'))
s0,f=R['W1_2']; kind,dat=f[0]
def ev(d,r0,r2,ta,sg): return mp.fsum(int(c)*r0**int(e[0])*r2**int(e[1])*ta**int(e[3])*sg**int(e[4]) for e,c in d.items())
for v in ('0.01','0.001','0.0001'):
    r0=mp.mpf('0.969'); r2=1/mp.mpf(v); ta=mp.mpf('0.22'); sg=mp.mpf('0.22')
    r3=-(r0+1+r2); rs=[r0,1,r2,r3]; E2=mp.fsum(rs[i]*rs[j] for i in range(4) for j in range(i+1,4)); E4=r0*r2*r3; D=E2*E2-4*E4
    F1,Y,P=[ev(dat[i],r0,r2,ta,sg) for i in range(3)]
    print(v,mp.nstr(F1,4),mp.nstr(Y,4),mp.nstr(P,4),'rel P',mp.nstr(P/(F1*F1*D),4),'val',mp.nstr((F1*mp.sqrt(D)-Y)/2,4))

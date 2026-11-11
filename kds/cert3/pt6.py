import pickle, mpmath as mp
mp.mp.dps=100
R=pickle.load(open('margin/zp_res.pkl','rb'))
s0,facs=R['Zp2']; kind,dat=facs[0]
def ev(d,r0,r2,ta): return mp.fsum(int(c)*r0**int(e[0])*r2**int(e[1])*ta**int(e[3]) for e,c in d.items())
for e_ in ('1e-3','1e-6'):
  e=mp.mpf(e_)
  for w in ('1e-2','1e-6'):
    r0=1-e; r2=1/(mp.mpf('0.486')*e); ta=1-e*mp.mpf(w)
    r3=-(r0+1+r2); rs=[r0,1,r2,r3]; E2=mp.fsum(rs[i]*rs[j] for i in range(4) for j in range(i+1,4)); E4=r0*r2*r3; D=E2*E2-4*E4
    F1,Y,P=[ev(dat[i],r0,r2,ta) for i in range(3)]
    raw=dat[6] if len(dat)>6 else None
    print(e_,w,'F1',mp.nstr(F1,4),'Y',mp.nstr(Y,4),'val',mp.nstr((F1*mp.sqrt(D)-Y)/2,4),'relP',mp.nstr(P/(F1*F1*D),4))

import pickle, mpmath as mp
mp.mp.dps=60
T=pickle.load(open('margin/scg_tau.pkl','rb'))
def ev(d,r0,r2): return mp.fsum(int(c)*mp.mpf(r0)**int(e[0])*mp.mpf(r2)**int(e[1]) for e,c in d.items())
for (s,be) in [(mp.mpf('0.176575'),mp.mpf('0.82755')),(mp.mpf('0.29920'),mp.mpf('0.05645'))]:
    r0=1-s/8; r2=1/(1-s*be/8); r3=-(r0+1+r2); rs=[r0,1,r2,r3]
    E2=mp.fsum(rs[i]*rs[j] for i in range(4) for j in range(i+1,4)); E4=r0*r2*r3; D=E2*E2-4*E4
    for c in ('q0','q1','v1'):
        F1,Y,P=[ev(T[('Q2',c)][i],r0,r2) for i in range(3)]
        print(float(s),c,'F1',mp.nstr(F1,5),'Y',mp.nstr(Y,5),'value',mp.nstr((F1*mp.sqrt(D)-Y)/2,5))

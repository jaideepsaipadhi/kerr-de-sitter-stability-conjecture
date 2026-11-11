exec(open("sal_build.py").read().split("\nres={}")[0])
import mpmath as mp
e=ex['Q2']
pts=[(mp.mpf('0.9858'),1/mp.mpf('0.8643'),mp.mpf('0.9763'))]
for r0v,r2v,tv in pts:
    r3v=-(r0v+1+r2v); rs=[r0v,1,r2v,r3v]
    E2=mp.fsum(rs[i]*rs[j] for i in range(4) for j in range(i+1,4)); E4=r0v*r2v*r3v
    bp=(E2+mp.sqrt(E2*E2-4*E4))/2
    for nm,P in (('n',e.n),('d',e.d)):
        c,fl=P.factor()
        for f,m in fl:
            val=mp.fsum(int(cc)*r0v**int(ee[0])*r2v**int(ee[1])*bp**int(ee[2])*tv**int(ee[3]) for ee,cc in f.to_dict().items())
            sc=mp.fsum(abs(int(cc))*abs(r0v)**int(ee[0])*abs(r2v)**int(ee[1])*abs(bp)**int(ee[2]) for ee,cc in f.to_dict().items())
            print(nm,f.degrees(),m,mp.nstr(val,4),'rel',mp.nstr(val/sc,3))

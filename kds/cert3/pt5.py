exec(open('sbe2_build.py').read().split("\nres={}")[0])
import mpmath as mp
mp.mp.dps=80
e=ex['W1_2']; cn,fl=e.n.factor()
for f,m in fl:
    if f.degrees()[2]==0: continue
    for v in ('0.01','0.0001'):
        r0=mp.mpf('0.969'); r2v=1/mp.mpf(v); ta=mp.mpf('0.22'); sgv=mp.mpf('0.22')
        r3=-(r0+1+r2v); rs=[r0,1,r2v,r3]; E2=mp.fsum(rs[i]*rs[j] for i in range(4) for j in range(i+1,4)); E4=r0*r2v*r3
        bp=2*E4/(E2-mp.sqrt(E2*E2-4*E4))
        terms={}
        for ee,c in f.to_dict().items():
            k=int(ee[2]); terms[k]=terms.get(k,0)+int(c)*r0**int(ee[0])*r2v**int(ee[1])*ta**int(ee[3])*sgv**int(ee[4])*bp**k
        tot=mp.fsum(terms.values()); mx=max(abs(x) for x in terms.values())
        print(f.degrees(),v,'f(b+)',mp.nstr(tot,4),'max term',mp.nstr(mx,4),'rel',mp.nstr(tot/mx,3))

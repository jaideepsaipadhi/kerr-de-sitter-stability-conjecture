exec(open("scg_build.py").read().split("\nres={}")[0])
import math
e=ex['H2']; c,fl=e.d.factor()
for (s,be) in [(0.066589,0.69766)]:
    r0v=1-s/8; vv=1-s*be/8; r2v=1/vv; r3v=-(r0v+1+r2v); rs=[r0v,1,r2v,r3v]
    E2=sum(rs[i]*rs[j] for i in range(4) for j in range(i+1,4)); E4=r0v*r2v*r3v
    bp=(E2+math.sqrt(E2*E2-4*E4))/2; bm=(E2-math.sqrt(E2*E2-4*E4))/2
    for f,m in fl:
        val=lambda bb: sum(float(int(cc))*r0v**int(ee[0])*r2v**int(ee[1])*bb**int(ee[2]) for ee,cc in f.to_dict().items())
        print(s,be,f.degrees(),'at b+',val(bp),'b-',val(bm))
    print('b+',bp,'b-',bm)

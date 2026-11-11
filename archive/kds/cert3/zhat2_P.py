import flint, pickle, sys, time
from fullcert import tocube
from ch2 import ctx4
from fastbern import pos_open_cube
k=sys.argv[1]
F1,Y,P=[ctx4.from_dict(x) for x in pickle.load(open('zhat2_res_%s.pkl'%k,'rb'))[k]]
t0=time.time(); g=F1.gcd(Y); print('gcd',g.degrees(),round(time.time()-t0,1),flush=True)
F1r=F1//g; Yr=Y//g
r0,r2,bb,tau=ctx4.gens(); r3=-(r0+1+r2)
e2=r0+r0*r2+r0*r3+r2+r3+r2*r3; e4=r0*r2*r3
Pr=F1r*F1r*(e2*e2-4*e4)-Yr*Yr
print('Pr',Pr.degrees(),len(Pr),round(time.time()-t0,1),flush=True)
q=tocube(Pr); print('cube',q.degrees(),round(time.time()-t0,1),flush=True)
print('pos',pos_open_cube(q),round(time.time()-t0,1),flush=True)

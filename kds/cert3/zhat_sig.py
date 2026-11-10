import flint, pickle, sys, numpy as np
from fractions import Fraction as Fr
from fullcert import tocube
from ch2 import ctx4
k,which=sys.argv[1],sys.argv[2]
F1,Y,P=pickle.load(open('zhat_res_%s.pkl'%k,'rb'))[k]
c,fl=ctx4.from_dict({'F1':F1,'P':P}[which]).factor()
big=max(fl,key=lambda fm:sum(fm[0].degrees()))[0]; q=tocube(big)
rng=np.random.default_rng(0); sg=[]; dd=q.to_dict()
for _ in range(60):
    D=10**6; x,v,t=[int(a*D) for a in rng.uniform(0.001,0.999,3)]
    dx,dv,dt=q.degrees()
    val=sum(int(cf)*x**int(e[0])*D**(dx-int(e[0]))*v**int(e[1])*D**(dv-int(e[1]))*t**int(e[2])*D**(dt-int(e[2])) for e,cf in dd.items())
    sg.append(1 if val>0 else -1)
print(k,which,'signs',sum(s>0 for s in sg),sum(s<0 for s in sg))

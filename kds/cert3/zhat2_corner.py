import flint, pickle, sys, itertools
from fractions import Fraction as Fr
from fullcert import tocube
from ch2 import ctx4
k=sys.argv[1]
F1,Y,_=[ctx4.from_dict(x) for x in pickle.load(open('zhat2_res_%s.pkl'%k,'rb'))[k]]
c,fl=F1.factor(); big=max(fl,key=lambda fm:sum(fm[0].degrees()))[0]; q=tocube(big); dd=q.to_dict()
def val(x,v,t): return sum(int(cf)*x**int(e[0])*v**int(e[1])*t**int(e[2]) for e,cf in dd.items())
for x in (Fr(1,2),Fr(3,4),Fr(9,10),Fr(99,100),Fr(1)):
    print(x,[ (v,t,'+' if val(x,v,t)>0 else ('0' if val(x,v,t)==0 else '-')) for v in (0,Fr(1,100),Fr(1,2),1) for t in (0,Fr(1,100),Fr(1,2),Fr(99,100),1)])

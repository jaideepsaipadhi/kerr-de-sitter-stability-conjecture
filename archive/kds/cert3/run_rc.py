import sys, time, pickle, flint
from ch2 import core3, res
from rcert import *
which,nm,idx=sys.argv[1],sys.argv[2],int(sys.argv[3])
F1,Y,P=res[which]; p3=core3(F1 if nm=='F1' else P)[idx]
p=Q.from_dict({e:flint.fmpq(int(c)) for e,c in p3.to_dict().items()})
import rcert; rcert.DEAD[0]=time.time()+50
t0=time.time(); ok=certify(p,tmax=45,label=f'{which}{nm}{idx}')
print(which,nm,idx,p3.degrees(),'CERTIFIED' if ok else 'FAILED',round(time.time()-t0,1)); print(LOG[-6:])

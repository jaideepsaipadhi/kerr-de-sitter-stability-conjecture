from ch2 import core3,res
from domcert import *
import sys,time
t0=time.time()
which,nm,idx=sys.argv[1],sys.argv[2],int(sys.argv[3]); eps=Fraction(1,int(sys.argv[4]))
F1,Y,P=res[which]; p=core3(F1 if nm=='F1' else P)[idx]
print(' E2', dominate_adaptive(p,(0,0,1),(1,2),0,(Fraction(0),1-eps),eps),round(time.time()-t0,1),flush=True)
print(' E1', dominate_adaptive(p,(1,0,1),(0,2),1,(eps,1-eps),eps),round(time.time()-t0,1),flush=True)
print(' E3', dominate_adaptive(p,(1,1,0),(0,1),2,(eps,Fraction(1)),eps),round(time.time()-t0,1),flush=True)

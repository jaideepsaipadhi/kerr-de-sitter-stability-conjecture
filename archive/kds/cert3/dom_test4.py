from ch2 import core3,res
from domcert import *
import sys
which,nm,idx=sys.argv[1],sys.argv[2],int(sys.argv[3])
F1,Y,P=res[which]; p=core3(F1 if nm=='F1' else P)[idx]
for e in (4,8,16):
  for cor,fl in (('x1v0t1',(1,0,1)),('x1v1t1',(1,1,1)),('x1v1t0',(1,1,0))):
    print(e,cor, dominate(p,fl,(0,1,2),None,None,Fraction(1,e))[0])
d=Fraction(1,8)
print('E3',dominate_adaptive(p,(1,1,0),(0,1),2,(Fraction(1,64),1-d),Fraction(1,64),maxd=10)[0])

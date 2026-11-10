from ch2 import core3,res
from domcert import *
import sys
which,nm,idx=sys.argv[1],sys.argv[2],int(sys.argv[3]); eps=Fraction(1,int(sys.argv[4])); d=Fraction(1,int(sys.argv[5]))
F1,Y,P=res[which]; p=core3(F1 if nm=='F1' else P)[idx]
print('E2',dominate_adaptive(p,(0,0,1),(1,2),0,(Fraction(0),1-d),eps,maxd=10)[0])
print('E1',dominate_adaptive(p,(1,0,1),(0,2),1,(d,1-d),eps,maxd=10)[0])
print('E3',dominate_adaptive(p,(1,1,0),(0,1),2,(d,Fraction(1)),eps,maxd=10)[0])

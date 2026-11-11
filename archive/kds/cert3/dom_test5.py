from ch2 import core3,res
from domcert import *
F1,Y,P=res['Q2']; p=core3(P)[0]
d=Fraction(1,16)
for e in (64,256,1024):
  r=dominate_adaptive(p,(1,0,1),(0,2),1,(d,1-d),Fraction(1,e),maxd=12)
  print('E1',e,r[0],r[1])
for e in (64,256):
  print('E3',e,dominate_adaptive(p,(1,1,0),(0,1),2,(Fraction(1,64),1-d),Fraction(1,e),maxd=10)[0])
  print('E2',e,dominate_adaptive(p,(0,0,1),(1,2),0,(Fraction(0),1-d),Fraction(1,e),maxd=10)[0])

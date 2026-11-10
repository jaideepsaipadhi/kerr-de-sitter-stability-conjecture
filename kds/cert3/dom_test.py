from ch2 import core3,res
from domcert import *
import sys
which,nm,idx=sys.argv[1],sys.argv[2],int(sys.argv[3])
F1,Y,P=res[which]; p=core3(F1 if nm=='F1' else P)[idx]
for eps in (Fraction(1,16),Fraction(1,64),Fraction(1,256)):
    print('eps',eps)
    # edge E2 {v=0,tau=1}: trans (v,T) ; u = x in [0, 1-eps]
    print(' E2', dominate(p,(0,0,1),(1,2),0,(0,1-eps),eps))
    # edge E1 {x=1,tau=1}: trans (X,T); u=v in [eps,1-eps]
    print(' E1', dominate(p,(1,0,1),(0,2),1,(eps,1-eps),eps))
    # edge E3 {x=1,v=1}: trans (X,V); u=tau in [eps,1-eps]? (tau free)
    print(' E3', dominate(p,(1,1,0),(0,1),2,(0,1-eps),eps))
    for cor,fl in (('x1v0t1',(1,0,1)),('x1v1t1',(1,1,1)),('x1v1t0',(1,1,0))):
        print(' corner',cor, dominate(p,fl,(0,1,2),None,None,eps))

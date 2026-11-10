import sys; from assemble import *; from remcert import run_resume
from ch2 import core3,res
w,nm,i=sys.argv[1],sys.argv[2],int(sys.argv[3]); cor=eval(sys.argv[4])
F1,Y,P=res[w]; p=core3(F1 if nm=='F1' else P)[i]
print(run_resume(p,regions(p,cor,f'{w}_{nm}_{i}'),f'st_{w}_{nm}_{i}.pkl',tmax=int(sys.argv[5])))

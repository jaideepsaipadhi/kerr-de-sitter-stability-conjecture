from ch2 import core3,res
from autocert import *
import sys,pickle
w,nm,i,k=sys.argv[1],sys.argv[2],int(sys.argv[3]),int(sys.argv[4]); cor=eval(sys.argv[5])
F1,Y,P=res[w]; p=core3(F1 if nm=='F1' else P)[i]
u,a,b,fa,fb=edges(p)[k]
cl=[0]*3;cl[a]=fa;cl[b]=fb; ch=list(cl); ch[u]=1
lo=F(1,cor.get(tuple(cl),10**9)); hi=F(1,cor.get(tuple(ch),10**9))
fl=[0]*3; fl[a]=fa; fl[b]=fb
for e in map(int,sys.argv[6].split(',')):
    r=dominate_adaptive(p,tuple(fl),(a,b),u,(lo,F(1)-hi),F(1,e),maxd=12)
    print((u,a,b,fa,fb),e,r[0],r[1],flush=True)
    if r[0]: break

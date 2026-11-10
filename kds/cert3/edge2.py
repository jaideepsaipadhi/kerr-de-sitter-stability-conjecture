from ch2 import core3,res
from autocert import *
from domcert import dom_tubes,dom_tubes2 as dom_tubes
import sys,pickle,time
w,nm,i,k=sys.argv[1],sys.argv[2],int(sys.argv[3]),int(sys.argv[4]); cor=eval(sys.argv[5])
F1,Y,P=res[w]; p=core3(F1 if nm=='F1' else P)[i]
u,a,b,fa,fb=edges(p)[k]
cl=[0]*3;cl[a]=fa;cl[b]=fb; ch=list(cl); ch[u]=1
lo=F(1,cor.get(tuple(cl),10**9)); hi=F(1,cor.get(tuple(ch),10**9))
fl=[0]*3; fl[a]=fa; fl[b]=fb
t=time.time()
r=dom_tubes(p,tuple(fl),(a,b),u,(lo,F(1)-hi),[F(1,int(x)) for x in sys.argv[6].split(",")],maxd=int(sys.argv[7]))
print((u,a,b,fa,fb),None if r is None else [(str(x[0][0]),str(x[0][1]),str(x[1])) for x in r],round(time.time()-t,1))
pickle.dump(r,open(f'tube_{w}_{nm}_{i}_{k}.pkl','wb'))

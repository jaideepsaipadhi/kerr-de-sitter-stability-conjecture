# assemble corner/edge regions from tube pickles and run remainder B&B
import sys,pickle,time
from fractions import Fraction as F
from autocert import edges
import remcert
def regions(p,cor,tag):
    regs=[]
    for c,e in cor.items(): regs.append(tuple((F(1)-F(1,e),F(1)) if f else (F(0),F(1,e)) for f in c))
    for k,(u,a,b,fa,fb) in enumerate(edges(p)):
        tubes=pickle.load(open(f'tube_{tag}_{k}.pkl','rb')); assert tubes is not None
        for (lo,hi),e in tubes:
            r=[None]*3; r[a]=(F(1)-e,F(1)) if fa else (F(0),e); r[b]=(F(1)-e,F(1)) if fb else (F(0),e); r[u]=(lo,hi); regs.append(tuple(r))
    return regs
if __name__=='__main__':
    from ch2 import core3,res
    w,nm,i=sys.argv[1],sys.argv[2],int(sys.argv[3]); cor=eval(sys.argv[4])
    F1,Y,P=res[w]; p=core3(F1 if nm=='F1' else P)[i]
    print(remcert.run(p,regions(p,cor,f'{w}_{nm}_{i}'),tmax=int(sys.argv[5])))

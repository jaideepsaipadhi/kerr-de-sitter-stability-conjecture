# Automatic certificate: p>0 on open cube. Detect vanishing cube edges, dominate corners/edges, B&B the rest.
import sys,time,itertools
from fractions import Fraction as F
from domcert import dominate,dominate_adaptive,c3,GX,GV,GT,ONE
from remcert import run
G=(GX,GV,GT); Z=0*GX
def edges(p):
    out=[]
    for u in range(3):
        a,b=[i for i in range(3) if i!=u]
        for fa,fb in itertools.product((0,1),repeat=2):
            sub=list(G); sub[a]=ONE if fa else Z; sub[b]=ONE if fb else Z
            if p.compose(*sub,ctx=c3).is_zero(): out.append((u,a,b,fa,fb))
    return out
def cert(p,verbose=True):
    E=edges(p); corners={}
    for u,a,b,fa,fb in E:
        for fu in (0,1):
            c=[0]*3; c[a]=fa;c[b]=fb;c[u]=fu; corners[tuple(c)]=None
    for c in corners:
        for e in (16,64,256,1024):
            if dominate(p,c,(0,1,2),None,None,F(1,e))[0]: corners[c]=F(1,e); break
        else: return dict(ok=False,corner=c)
    regs=[]
    for c,e in corners.items():
        regs.append(tuple((F(1)-e,F(1)) if f else (F(0),e) for f in c))
    for u,a,b,fa,fb in E:
        cl=[0]*3;cl[a]=fa;cl[b]=fb;cl[u]=0; ch=list(cl); ch[u]=1
        lo=corners.get(tuple(cl)) or F(0); hi=corners.get(tuple(ch)) or F(0)
        fl=[0]*3; fl[a]=fa; fl[b]=fb
        for e in (64,256,1024,4096):
            if dominate_adaptive(p,tuple(fl),(a,b),u,(lo,F(1)-hi),F(1,e),maxd=12)[0]: break
        else: return dict(ok=False,edge=(u,a,b,fa,fb))
        r=[None]*3; r[a]=(F(1)-e,F(1)) if fa else (F(0),e); r[b]=(F(1)-e,F(1)) if fb else (F(0),e); r[u]=(lo,F(1)-hi)
        regs.append(tuple(r))
    if verbose: print(' edges',E,' corners',{k:str(v) for k,v in corners.items()},flush=True)
    return run(p,regs)
if __name__=='__main__':
    from ch2 import core3,res
    w,nm=sys.argv[1],sys.argv[2]
    F1,Y,P=res[w]
    for i,p in enumerate(core3(F1 if nm=='F1' else P)):
        t=time.time(); print(w,nm,i,p.degrees(),cert(p),round(time.time()-t,1),flush=True)

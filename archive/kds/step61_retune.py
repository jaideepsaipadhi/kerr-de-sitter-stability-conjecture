import sys
from step32_four_param import check
def grid(delta, beta=2.0):
    worst=[-9,-9,None,None]; fails=0; n=0
    for t in (1.01,1.1,1.5,2,5,20,100,1e3,1e4,1e5):
      for zf in (0.05,0.2,0.5,0.9,0.95):
        x0=zf*(t-1); z0=1+x0; Bf=5.5*t-1.25
        for Bm in (1,10,1e3):
          for ratio in (0.01,0.1,1,10,100,1e3,1e4,1e6):
            B=Bm*Bf; A=ratio*B; L1=A+B-3
            ok,info=check(t,z0,A,B,delta,0.5*min(0.9,2*x0*L1/(t-1)),beta,zdecay=True)
            if isinstance(info,str): continue
            n+=1; fails+= (not ok)
            if info[0]>worst[0]: worst[0]=info[0]; worst[2]=(t,zf,Bm,ratio)
            if info[1]>worst[1]: worst[1]=info[1]; worst[3]=(t,zf,Bm,ratio)
    return n,fails,worst
d=float(sys.argv[1]); b=float(sys.argv[2]) if len(sys.argv)>2 else 2.0
n,f,w=grid(d,b); print(f'delta={d} beta={b}: {n} cases, fails {f}, worst C {w[0]:+.4f} at {w[2]}, worst D {w[1]:+.4f} at {w[3]}', flush=True)

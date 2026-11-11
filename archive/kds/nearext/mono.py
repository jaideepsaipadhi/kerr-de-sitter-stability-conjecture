import numpy as np
exec(open('/tmp/lz0c.py').read().split("def F(")[0])
def Lz0(g,fr,m):
    Q2,Q1,Q0,z0,A2,C0,t=coeffs(g,fr); return m*m*Q2+m*Q1+Q0+(4*C0 if m==1 else 0)
nm=0; nc=0; tot=0; worst=9
for rho in np.concatenate([np.linspace(0.02,0.98,30),1-np.logspace(-2,-6,8)]):
  for v in np.concatenate([np.logspace(-5,-1,8),np.linspace(0.1,0.995,30)]):
    try: g=geom(rho,v)
    except Exception: continue
    fr=np.linspace(0,1,201)
    for m in (1,2,3,6):
        y=np.array([Lz0(g,f,m) for f in fr]); tot+=1
        if np.any(np.diff(y)<-1e-12*abs(y).max()): nm+=1
        if np.any(np.diff(y,2)>1e-9*abs(y).max()): nc+=1
        worst=min(worst,y[0]/g['C0'])
print('tot',tot,'nonmonotone',nm,'nonconcave',nc,'min L(z0)/C0 at w2',worst)

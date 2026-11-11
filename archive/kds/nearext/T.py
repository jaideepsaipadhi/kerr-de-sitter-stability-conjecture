import numpy as np
exec(open('/tmp/lz0c.py').read().split("def coeffs")[0])
def T(g,w,m):
    P,Om,a,al,Xi,C0=g['P'],g['Om'],g['a'],g['al'],g['Xi'],g['C0']
    u0=-P[0]*(w-Om[0]); u3=P[3]*(w-Om[3])
    d=1 if m==1 else 0
    return 4*Xi**2*u0*u3/C0+1+1/m-3/m**2+4*d/m**2-5*al/m**2-4*al*Xi/m-2*Xi**2*a*w*(1+2/m)
worst=(9,)
for rho in np.concatenate([np.linspace(0.01,0.99,40),1-np.logspace(-2,-7,12)]):
  for v in np.concatenate([np.logspace(-6,-1,12),np.linspace(0.1,0.999,40)]):
    try: g=geom(rho,v)
    except Exception: continue
    for w in (g['w1'],g['w2']):
      for m in (1,2,3,5,10,100,10**6):
        val=T(g,w,m)
        if val<worst[0]: worst=(val,rho,v,'w1' if w==g['w1'] else 'w2',m)
print(worst)

# Normal hyperbolicity of KN trapping: at every trapped orbit (R=R'=0), check R''>0 (radial instability).
import numpy as np, random
from scipy.optimize import brentq
bad=0; worst=1e9
for _ in range(3000):
    a=random.uniform(-.999,.999); q=random.uniform(0,(1-a*a)**.5*.999); M=1
    rp=1+(1-a*a-q*q)**.5
    for W in np.concatenate([np.linspace(-3,3,200)]):   # Omega=E/L ; L=1 (also L->0 limit below)
        D=lambda r:r*r-2*r+a*a+q*q
        P=lambda r:(r*r+a*a)*W-a
        G=lambda r:4*r*W*D(r)-2*(r-1)*P(r)
        xs=np.linspace(rp*1.0001,40,2000); g=G(xs)
        for i in np.where(np.sign(g[:-1])!=np.sign(g[1:]))[0]:
            r=brentq(G,xs[i],xs[i+1]); K=P(r)**2/D(r)
            if K<(1-a*W)**2-1e-12: continue        # not a real geodesic
            R2=2*(2*r*W)**2+2*P(r)*2*W-2*K       # R''=2P'^2+2P P''-D''K
            worst=min(worst,R2/W**2 if W else R2); bad+= R2<=0
print("orbits with R''<=0:",bad,"  min normalized R'':",worst)

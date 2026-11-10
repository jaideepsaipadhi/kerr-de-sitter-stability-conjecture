import numpy as np, mpmath as mp, sympy as sp
from scipy.integrate import solve_ivp
exec(open('step138.py').read().split("mp.mp.dps")[0].replace("import sympy as sp, numpy as np, mpmath as mp, sys",""))
mp.mp.dps=30
for tt,zf,AA,bb in ((57.565,0.053,0.8037341928942603,0.49849745927989964),(430.069,0.016,16.06666690380505,1.8310315893064912)):
    zz=1+zf*(tt-1); x0=zz-1; BB=3*zz-AA/zz+bb
    c=[mp.mpf(g(tt,zz,AA,BB)) for g in cf]
    r=[mp.re(x) for x in mp.polyroots(c,maxsteps=300,extraprec=300) if abs(mp.im(x))<1e-15 and 0<mp.re(x)<x0]
    Sh=float(r[0]); zh=zz-Sh; mu_=zh*(zh-1)*(tt-zh); mup=(zh-1)*(tt-zh)+zh*(tt-zh)-zh*(zh-1); H=2*mu_+Sh*mup; w0=H/Sh
    Lf=lambda q: AA/q+BB-3*q; muf=lambda q:q*(q-1)*(tt-q)
    sol=solve_ivp(lambda q,w:[(-w[0]**2+2*(zz-q)*Lf(q)*w[0]-2*muf(q)*(zz-q)*(AA/q**2+3))/muf(q)],(zh,zz),[w0],rtol=1e-11,atol=1e-14,dense_output=True)
    qq=np.linspace(zh,zz,9); print('w/w0',np.round(sol.sol(qq)[0]/w0,3),'2sL/w0',np.round(2*(zz-qq)*Lf(qq)/w0,3), 'mu0',muf(zz),'Sh',Sh,'x0',x0)
    for fac in (1,1.5,3,10):
        M=fac*w0
        sv=solve_ivp(lambda q,v:[((2*(zz-q)*Lf(q)-M)*v[0])/muf(q)-2*(zz-q)*(AA/q**2+3)],(zh,zz),[w0],rtol=1e-11,atol=1e-14,dense_output=True)
        vv=sv.sol(np.linspace(zh,zz,200))[0]; print(' fac',fac,'vmax/M',round(vv.max()/M,3),'v_end/w0',round(vv[-1]/w0,4))

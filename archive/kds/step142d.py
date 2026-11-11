exec(open('step142c.py').read().split("best={}")[0])
from scipy.integrate import solve_ivp
n=0
for P in out:
    if condmax(P,1.5,0.5): continue
    tt,zz,AA,BB,Sh=P; zh=zz-Sh; mu_=zh*(zh-1)*(tt-zh); mup=(zh-1)*(tt-zh)+zh*(tt-zh)-zh*(zh-1); w0=(2*mu_+Sh*mup)/Sh
    Lf=lambda q:AA/q+BB-3*q; muf=lambda q:q*(q-1)*(tt-q)
    sol=solve_ivp(lambda q,w:[(-w[0]**2+2*(zz-q)*Lf(q)*w[0]-2*muf(q)*(zz-q)*(AA/q**2+3))/muf(q)],(zh,zz),[w0],method='LSODA',rtol=1e-10,atol=1e-12,dense_output=True)
    sg=np.array([1,.99,.97,.95,.9,.8,.6,.4,.2,0.05]); q=zz-Sh*sg
    print('t%.0f x0 %.3f xh %.2e'%(tt,zz-1,zh-1),'w/w0',np.round(sol.sol(q)[0]/w0,3),' 2sL/w0',np.round(2*(zz-q)*Lf(q)/w0,3))
    n+=1
    if n>=3: break

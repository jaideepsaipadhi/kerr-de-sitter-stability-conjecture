exec(open('sbe_build.py').read().split("\nres={}")[0])
G_=Gam; rho2=F(21,10*one)*(-Hm[2])*lamv/den
import numpy as np, warnings, os; warnings.filterwarnings('ignore')
os.chdir('../nearext'); exec(open('cos176.py').read().split("rng=np.random.default_rng(2)")[0])
for (r0v,r2v,fr) in [(0.5,3.0,0.3),(0.9,1.5,0.7)]:
    Lm,av,alv,w1=kds(r0v,1.0,r2v)
    t_,z0_,A_,B_,q1_=pdq(r0v,1.0,r2v,2,fr); Lt=A_/t_+B_-3*t_; T_=t_-z0_
    K_=1.05*(-Lt)*t_*(t_-1)/(2*T_*T_)
    pe=lambda p: sum(float(int(c))*r0v**int(e[0])*r2v**int(e[1])*(av*av)**int(e[2])*fr**int(e[3]) for e,c in p.to_dict().items())
    ev=lambda Fo: pe(Fo.n)/pe(Fo.d)
    print('Gamma',ev(G_),A_*T_/q1_**2,' rho2',ev(rho2),K_*T_/q1_**2)

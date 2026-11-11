import flint
exec(open('sal_build.py').read().split("exec(\"def bcoeffs")[0])
HH=ex['H2']
import numpy as np, warnings, os; warnings.filterwarnings('ignore')
os.chdir('../nearext'); exec(open('cos176.py').read().split("rng=np.random.default_rng(2)")[0])
exec(open('../cert3/imchk.py').read().split("FF=(f1,f2,Xi,lamv)")[1].split("os.chdir")[1].split("\n",1)[1].split("for (r0v")[0])
for (r0v,r2v,fr) in [(0.5,3.0,0.3),(0.9,1.5,0.7),(0.2,20.0,0.5)]:
    Lm,av,alv,w1=kds(r0v,1.0,r2v); lt=Lm/3
    t_,z0_,A_,B_,q1_=pdq(r0v,1.0,r2v,2,fr); Lt=A_/t_+B_-3*t_; u,et,_=ex(r0v,1.0,r2v,2,fr)
    F=-(u+t_*et); W=3*t_*(t_-1)*et*et*F+2.1*Lt*(2*F+(t_-1)*(F+3*(-u)))
    pe=lambda p: sum(float(int(c))*r0v**int(e[0])*r2v**int(e[1])*(av*av)**int(e[2])*fr**int(e[3]) for e,c in p.to_dict().items())
    print(pe(HH.n)/pe(HH.d), lt*W/(2*av)/4)

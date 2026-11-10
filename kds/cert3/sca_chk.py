import flint
exec(open('sca_build.py').read().split("exec(\"def bcoeffs")[0])
EXC=ex['CH2']; EXA=ex['AH2']
import numpy as np, warnings, os; warnings.filterwarnings('ignore')
os.chdir('../nearext'); exec(open('cos176.py').read().split("rng=np.random.default_rng(2)")[0])
for (r0v,r2v,fr) in [(0.5,3.0,0.3),(0.9,1.5,0.7),(0.2,20.0,0.5)]:
    Lm,av,alv,w1=kds(r0v,1.0,r2v); lt=Lm/3
    t_,z0_,A_,B_,q1_=pdq(r0v,1.0,r2v,2,fr); Lt=A_/t_+B_-3*t_; T_=t_-z0_; a_=z0_-1
    pe=lambda p: sum(float(int(c))*r0v**int(e[0])*r2v**int(e[1])*(av*av)**int(e[2])*fr**int(e[3]) for e,c in p.to_dict().items())
    ev=lambda Fo: pe(Fo.n)/pe(Fo.d)
    dC=(2*lt*q1_**2*T_-1.05*t_*(t_-1)*(-lt*Lt))/4
    dA=(6*lt*q1_**2*T_**2/(2+4*a_+T_)-1.05*t_*(t_-1)*(-lt*Lt))/4
    print(ev(EXC),dC,'|',ev(EXA),dA)

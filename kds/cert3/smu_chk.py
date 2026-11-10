import flint
exec(open('smu_build.py').read().split("exec(\"def bcoeffs")[0])
Z0F=z0; HH=ex["H2"]
import numpy as np, warnings; warnings.filterwarnings('ignore')
import os; os.chdir('../nearext'); exec(open('cos176.py').read().split("rng=np.random.default_rng(2)")[0])
for (r0v,r2v,fr) in [(0.5,3.0,0.3),(0.9,1.5,0.7),(0.2,20.0,0.5)]:
    Lm,av,alv,w1=kds(r0v,1.0,r2v)
    t_,z0_,A_,B_,q1_=pdq(r0v,1.0,r2v,2,fr)
    q=lambda f: flint.fmpq(int(round(f*10**12)),10**12)
    pe=lambda p: sum(float(int(c))*r0v**int(e[0])*r2v**int(e[1])*(av*av)**int(e[2])*fr**int(e[3]) for e,c in p.to_dict().items())
    def ev(Fo): return pe(Fo.n)/pe(Fo.d)
    print('H2',ev(HH),' direct',(Lm/3)*(A_/t_+B_-3*t_)/4, ' z0',ev(Z0F),z0_)

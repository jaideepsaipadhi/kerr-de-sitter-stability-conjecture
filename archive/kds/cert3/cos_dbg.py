import flint, numpy as np, mpmath as mp
exec(open('zhat_build.py').read().split("Z=t/R")[0])
mp.mp.dps=40
def ev(e,vals):
    def p(q):
        s=mp.mpf(0)
        for ex,c in q.to_dict().items(): s+=int(c)*vals[0]**int(ex[0])*vals[1]**int(ex[1])*vals[2]**int(ex[2])*vals[3]**int(ex[3])
        return s
    return p(e.n)/p(e.d)
f1=U[1]+U[3]; f2=U[2]+U[3]
import sys; sys.path.insert(0,'/home/claude/kerr-stability/kds/nearext')
r0=mp.mpf('0.9'); r2=mp.mpf('10'); fr=mp.mpf('0.5')
r3=-(r0+1+r2); e2v=r0+r0*r2+r0*r3+r2+r3+r2*r3; e4v=r0*r2*r3; bv=(e2v+mp.sqrt(e2v**2-4*e4v))/2
vals=(r0,r2,bv,fr)
print('f1',ev(f1,vals),'f2',ev(f2,vals),'U0',ev(U[0],vals),'U3',ev(U[3],vals),'Xi',ev(Xi,vals),'lamv',ev(lamv,vals),'t',ev(t,vals),'b',bv)

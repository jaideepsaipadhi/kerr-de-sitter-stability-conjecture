exec(open('step143.py').read().split("def two(P,th):")[0])
from scipy.integrate import trapezoid
def tail(P):
    tt,zz,AA,BB,Sh=P
    sg=np.concatenate([np.geomspace(1e-9,1,6000)]); q=zz-Sh*sg; s=Sh*sg; L=AA/q+BB-3*q; Lp=AA/q**2+3; mu=q*(q-1)*(tt-q)
    D=s*s*L*L-2*mu*s*Lp
    neg=np.where(D<=0)[0]
    if len(neg)==0: return True
    i=neg.max()           # tail = sg[0..i], last index where D<=0; switch at i (D~0)
    ws=s[i]*L[i]
    zq=q[:i+1][::-1]; f1=(1/mu[:i+1])[::-1]; f2=(2*s[:i+1]*Lp[:i+1]-s[:i+1]**2*L[:i+1]**2/mu[:i+1])[::-1]
    I1=trapezoid(f1,zq); I2=trapezoid(f2,zq)
    return ws-ws*ws*I1-I2>0
G=gen(1500,5)
print('corner',sum(tail(Q) for Q in out),len(out),'generic',sum(tail(Q) for Q in G),len(G))

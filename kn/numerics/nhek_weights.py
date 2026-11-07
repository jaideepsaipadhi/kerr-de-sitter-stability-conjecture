# Extremal Kerr (M=1,a=1), s=-2, omega = m*Omega_H^ext = m/2 (c = a*omega = m/2).
# Spin-weighted spheroidal eigenvalue A via spectral method; NHEK weight delta^2 = 7m^2/4 - (s+1/2)^2 - A.
# delta^2 > 0  <=> oscillatory weights (mu < -1/4): the near-extremal "balance" regime.
import numpy as np
from sympy.physics.wigner import clebsch_gordan as CG
def cosmat(s,m,L):
    ls=list(range(max(abs(s),abs(m)),L))
    n=len(ls); C=np.zeros((n,n))
    for i,l in enumerate(ls):
        for j,lp in enumerate(ls):
            if abs(l-lp)<=1:
                # <s l m|cos|s lp m> = sqrt((2lp+1)/(2l+1)) <lp m 1 0|l m><lp -s 1 0|l -s>
                C[i,j]=float(np.sqrt((2*lp+1)/(2*l+1))*CG(lp,1,l,m,0,m)*CG(lp,1,l,-s,0,-s))
    return ls,C
s=-2
for m in range(1,9):
    ls,C=cosmat(s,m,m+40)
    c=m/2
    H=np.diag([l*(l+1)-s*(s+1) for l in ls]) - c*c*(C@C) + 2*c*s*C
    A=np.sort(np.linalg.eigvals(H).real)
    d2=7*m*m/4-(s+0.5)**2-A[:4]
    print(f"m={m}:", "  ".join(f"l={ls[0]+k}: d2={d2[k]:+.3f}" for k in range(4)))

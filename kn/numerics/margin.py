# Lemma A margin: min_{r>r+} G(r)/(normaliser) as Omega -> omega_H and toward extremality.
import numpy as np
def margin(a,q,eta):
    rp=1+np.sqrt(1-a*a-q*q); wH=a/(rp*rp+a*a); W=wH*(1-eta)
    r=rp+np.logspace(-8,2,20000)
    D=r*r-2*r+a*a+q*q; P=(r*r+a*a)*W-a
    G=4*r*W*D-2*(r-1)*P
    return (G/(r**3*W+a)).min()          # scale-free normalisation
for a,q in [(.5,.3),(.8,.3),(.9,.4),(.95,.3)]:
    rp=1+np.sqrt(1-a*a-q*q)
    print(f"a={a} Q={q} r+-M={rp-1:.3f}: ",
          [f"eta={e}: {margin(a,q,e)/(e*(rp-1)):.3f}" for e in (1e-1,1e-2,1e-3,1e-4)])

import numpy as np
from scipy.special import eval_jacobi, roots_jacobi
from numpy.polynomial.legendre import leggauss
def lam_min(m,s,al,xi,nb=14,nq=400):
    a,b=abs(m+s),abs(m-s); l0=max(abs(m),abs(s))
    x,w=leggauss(nq)
    def phi(n,x):
        return (1-x)**(a/2)*(1+x)**(b/2)*eval_jacobi(n,a,b,x)
    h=1e-6
    Ph=np.array([phi(n,x) for n in range(nb)]); dP=np.array([(phi(n,x+h)-phi(n,x-h))/(2*h) for n in range(nb)])
    # normalize
    nrm=np.sqrt((Ph**2*w).sum(1)); Ph/=nrm[:,None]; dP/=nrm[:,None]
    VH=-(1+al)**2*xi**2/al+2*al*x**2-(1+al)/(1+al*x**2)*(2*s*(al*m-(1+al)*xi)*x-(1+al)**2*xi**2/al+2*m*(1+al)*xi+s**2)+(1+al)**2*m**2/((1+al*x**2)*(1-x**2))+(1+al)*(s**2+2*s*m*x)/(1-x**2)
    V0=(m+s*x)**2/(1-x**2)
    H=np.diag([(l0+n)*(l0+n+1)-s**2 for n in range(nb)]).astype(float)
    H+= (dP*(al*x**2*(1-x**2))*w)@dP.T + (Ph*((VH-V0))*w)@Ph.T
    return np.linalg.eigvalsh((H+H.T)/2)[0]
for m in ():
  for s in (-2,2):
    worst=min((lam_min(m,s,al,xi),al,xi) for al in np.linspace(1e-4,0.0718,6) for xi in np.linspace(0,m*0.5,11))
    print(m,s,worst)
print('check', lam_min(2,-2,1e-4,0.0), lam_min(2,2,1e-4,0.0), lam_min(1,-2,1e-4,0.0))
for xi in [0.1,0.2,0.3,0.5,0.7,1.0]: print(xi, lam_min(2,-2,1e-4,xi), lam_min(1,-2,1e-4,xi/2))

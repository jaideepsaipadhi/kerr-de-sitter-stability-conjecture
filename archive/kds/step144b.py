exec(open('step143.py').read().split("def two(P,th):")[0])
from scipy.integrate import cumulative_trapezoid as ct
def tail(P,th=0.95):
    tt,zz,AA,BB,Sh=P; zh=zz-Sh; mu_=zh*(zh-1)*(tt-zh); mup=(zh-1)*(tt-zh)+zh*(tt-zh)-zh*(zh-1); w0=(2*mu_+Sh*mup)/Sh
    sg=np.geomspace(1e-9,1,4000); q=zz-Sh*sg; s=Sh*sg; L=AA/q+BB-3*q; Lp=AA/q**2+3; mu=q*(q-1)*(tt-q)
    D=s*s*L*L-2*mu*s*Lp
    K=np.minimum(w0, s*L+th*np.sqrt(np.maximum(D,0)))
    # zone1 validity from start (sg=1) down to switch: D>0 and K nonincreasing in z (i.e. nondecreasing in sg) and K<kp
    kp=s*L+np.sqrt(np.maximum(D,0))
    okz=(D>0)&(K<kp+1e-15)
    # cumulative integrals from z0 (sg small) upward: z increases as sg decreases
    z_=q[::-1]  # increasing in z? q decreasing in sg -> reversed increasing
    # integral from z_s to z0 of g dz = integral over sg in (0, sg_s) of g*Sh dsg
    I1=ct(1/mu*Sh,sg,initial=0); I2=ct(np.maximum(-D,0)/mu*Sh,sg,initial=0)
    mono=np.concatenate([[True],np.diff(K)>=-1e-12*np.abs(K[1:])])
    # zone1 must be valid on [i, end]
    valid_from=np.ones(len(sg),bool)
    good=okz&mono
    suffix=np.flip(np.cumprod(np.flip(good)))>0
    for i in np.where(suffix)[0][::5]:
        vs=K[i]
        if vs-vs*vs*I1[i]-I2[i]>0: return True
    return False
G=gen(1500,5)
for th in (0.9,0.99):
    print(th,'corner',sum(tail(Q,th) for Q in out),len(out),'generic',sum(tail(Q,th) for Q in G),len(G),flush=True)

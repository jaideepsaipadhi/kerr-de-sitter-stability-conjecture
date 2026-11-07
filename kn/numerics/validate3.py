import sys; sys.argv=['x']
exec(open('realmode_scan.py').read().split("if __name__")[0])
av,Qv,mv=0.001,0.0001,2; N=8; full=setup(av,Qv,mv,N=N)
def mats(R,W):
    A=full(R,W); return [X[:N,:N] for X in A]
def smin1(W,x0=1e-5,rbig=60.,rmid=3.):
    rp=1+np.sqrt(1-av**2-Qv**2); rm_=1-np.sqrt(1-av**2-Qv**2); OH=av/(rp**2+av**2)
    sig=(W-mv*OH)*(rp**2+av**2)/(rp-rm_)
    def F(R,y):
        Z=y.view(complex).reshape(2,N,N); A0,A1,A2=mats(R,W)
        return np.stack([Z[1],-np.linalg.solve(A2,A1@Z[1]+A0@Z[0])]).ravel().view(float)
    q=np.full(N,2-1j*sig)
    Zh=np.stack([np.diag(x0**q),np.diag(q*x0**(q-1))])
    a1=solve_ivp(F,(rp+x0,rmid),Zh.ravel().view(float),rtol=1e-9,atol=1e-30).y[:,-1].view(complex).reshape(2,N,N)
    rs=lambda R:R+(rp**2+av**2)/(rp-rm_)*np.log(R-rp)-(rm_**2+av**2)/(rp-rm_)*np.log(R-rm_)
    ph=np.exp(1j*W*rs(rbig))*rbig**3*np.ones(N); dph=ph*(3/rbig+1j*W*(rbig**2+av**2)/((rbig-rp)*(rbig-rm_)))
    a2=solve_ivp(F,(rbig,rmid),np.stack([np.diag(ph),np.diag(dph)]).ravel().view(float),rtol=1e-9,atol=1e-30).y[:,-1].view(complex).reshape(2,N,N)
    Mm=np.block([[a1[0],a2[0]],[a1[1],a2[1]]]); Mm=Mm/np.linalg.norm(Mm,axis=0)
    return np.linalg.svd(Mm,compute_uv=False)[-3:]
for W in [0.37367168-0.08896232j,0.40-0.089j,0.35-0.089j,0.3737-0.06j,0.3737-0.12j]:
    print(W, smin1(W),flush=True)

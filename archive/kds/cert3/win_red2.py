import flint, time, pickle, sympy as sp, numpy as np
t0=time.time()
exec(open('win_flint.py').read().split("print('R z0'")[0])
exec("Xi=1-F(b*b,e4); al=F(-b*b,e4); lamv=F(-b,e4); K=F(one,(r1-r0)*(r2-r3))")
exec(open('win_flint.py').read().split("K=F(one,(r1-r0)*(r2-r3))")[1].split("out={}")[0].replace("print(","(lambda *a,**k:None)("))
# reduce N*D mod g(b) = b^2 - e2 b + e4  (monic in b): treat as poly in b
S=sp.symbols('r0 r2 b tau'); sr0,sr2,sb,stau=S
def tosp(p): return sp.sympify(str(p).replace('^','**'),locals=dict(zip(['r0','r2','b','tau'],S)))
E2,E4=tosp(e2),tosp(e4); G=sp.Poly(sb**2-E2*sb+E4,sb)
res={}
for k,e in ex.items():
    NP=tosp(e.n*e.d)
    rem=sp.Poly(NP,sb).rem(G); cs=rem.all_coeffs(); cs=[0]*(2-len(cs))+cs
    F1c,F0c=[sp.expand(c) for c in cs]
    Y=sp.expand(-(2*F0c+F1c*E2)); D=sp.expand(E2**2-4*E4); P=sp.expand(F1c**2*D-Y**2)
    res[k]=(F1c,Y,P); print(k,'reduced',round(time.time()-t0,1),flush=True)
pickle.dump(res,open('win_res.pkl','wb'))
# numeric sign survey on the physical domain
for k,(F1c,Y,P) in res.items():
    fF1=sp.lambdify((sr0,sr2,stau),F1c); fY=sp.lambdify((sr0,sr2,stau),Y); fP=sp.lambdify((sr0,sr2,stau),P)
    g0=np.linspace(0.02,0.98,25); vv=np.linspace(0.02,0.98,25); tt=np.linspace(0,1,11)
    A,B_,C=np.meshgrid(g0,1/vv,tt)
    s1=np.sign(fF1(A,B_,C)); sp_=np.sign(fP(A,B_,C)); sy=np.sign(fY(A,B_,C))
    print(k,'F1 signs',s1.min(),s1.max(),' P signs',sp_.min(),sp_.max(),' Y signs',sy.min(),sy.max())

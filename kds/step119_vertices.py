import pickle, sympy as sp, mpmath as mp, random
mp.mp.dps=60
t_,z0_,A_,B_=sp.symbols('t z0 A B',positive=True)
Tj=pickle.load(open('/home/claude/kds/T118.pkl','rb'))
f={j:sp.lambdify((t_,z0_,A_,B_),Tj[j],'mpmath') for j in (2,3,4,5,6)}
random.seed(3); worst={j:-mp.inf for j in (4,5,6)}; arg={}
cnt={'V2':0,'n':0}; c_ok=True; cbad=0
for _ in range(20000):
    t=1+mp.mpf(10)**random.uniform(-2,5); z0=1+mp.mpf(random.uniform(1e-4,1-1e-4))*(t-1); x0=z0-1
    Bf=mp.mpf(11)/2*t-mp.mpf(5)/4; Ac2=max(mp.mpf(0),z0*((t-1)/(2*x0*x0)-3)); P=3+(t-1)/(2*x0)
    V=[(max(Ac2,P-Bf),Bf)]
    if P-Ac2>Bf: V.append((Ac2,P-Ac2)); cnt['V2']+=1
    cnt['n']+=1
    sc=abs(f[6](t,z0,V[0][0],V[0][1]))
    for j in (4,5,6):
        v=max(f[j](t,z0,a,b) for a,b in V)/sc
        if v>worst[j]: worst[j]=v; arg[j]=(float(t),float((z0-1)/(t-1)))
print(cnt); print({j:(float(worst[j]),arg[j]) for j in worst})

import flint, pickle, sys, time
from fullcert import sign_on_cube, tocube
from ch2 import ctx4
k=sys.argv[1]; which=sys.argv[2]
F1,Y,P=pickle.load(open('zhat_res_%s.pkl'%k,'rb'))[k]
d={'F1':F1,'P':P}[which]
t0=time.time(); c,fl=ctx4.from_dict(d).factor(); print('factored',round(time.time()-t0,1),[(f.degrees(),m) for f,m in fl],flush=True)
s=1 if c>0 else -1
for f,m in fl:
    sg=sign_on_cube(f); print(f.degrees(),m,sg,round(time.time()-t0,1),flush=True)
    if sg is None: s=None; break
    if m%2: s*=sg
print(k,which,'=> sign',s)

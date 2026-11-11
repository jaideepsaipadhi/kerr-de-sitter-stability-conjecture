import flint, pickle, sys, time
from fullcert import tocube
from ch2 import ctx4
from fastbern import pos_open_cube
k=sys.argv[1]; which=sys.argv[2]
F1,Y,P=pickle.load(open('zhat2_res_%s.pkl'%k,'rb'))[k]
d={'F1':F1,'P':P,'Y':Y}[which]
t0=time.time(); c,fl=ctx4.from_dict(d).factor(); print('factored',round(time.time()-t0,1),flush=True)
s=1 if c>0 else -1
for f,m in fl:
    q=tocube(f)
    sg=1 if pos_open_cube(q) else (-1 if pos_open_cube(-q) else None)
    print(q.degrees(),m,sg,round(time.time()-t0,1),flush=True)
    if sg is None: s=None
    elif m%2: s*=sg
print(k,which,'=> sign',s)

import flint, pickle, sys
from fullcert import tocube
from ch2 import ctx4
from fastbern import pos_open_cube
k=sys.argv[1]
F1,Y,P=[ctx4.from_dict(x) for x in pickle.load(open('zhat2_res_%s.pkl'%k,'rb'))[k]]
g=F1.gcd(Y); c,fl=g.factor()
for f,m in fl:
    q=tocube(f); print(q.degrees(),m,1 if pos_open_cube(q) else (-1 if pos_open_cube(-q) else None))

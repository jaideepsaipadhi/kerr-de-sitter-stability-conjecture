import flint, pickle, sys, time
from fullcert import tocube
from ch2 import ctx4
from fastbern import subdiv_pos
k=sys.argv[1]
F1,Y,_=[ctx4.from_dict(x) for x in pickle.load(open('zhat2_res_%s.pkl'%k,'rb'))[k]]
c,fl=F1.factor(); big=max(fl,key=lambda fm:sum(fm[0].degrees()))[0]; q=tocube(big)
t0=time.time(); print(subdiv_pos(q,maxdepth=9,tmax=45), subdiv_pos(-q,maxdepth=2,tmax=5)[0], round(time.time()-t0,1))

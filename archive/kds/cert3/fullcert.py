# Sign certificate of F1 and P (for Q2,H2,H1) on the open window cube via factorisation + nonneg-Bernstein.
import sys,time
from ch2 import res,ctx4,c3
from nncert import run
def tocube(f):
    dd=f.to_dict(); assert all(e[2]==0 for e in dd), 'b present'
    n2=max(e[1] for e in dd)
    return c3.from_dict({(e[0],n2-e[1],e[3]):int(c) for e,c in dd.items()})
def sign_on_cube(f):
    q=tocube(f)
    r=run(q,None,tmax=50)
    if r['ok']: return 1
    r=run(-q,None,tmax=50)
    if r['ok']: return -1
    return None
for w in sys.argv[1:]:
    F1,Y,P=res[w]
    for nm,d in (('F1',F1),('P',P)):
        c,fl=ctx4.from_dict(d).factor(); s=1 if c>0 else -1; info=[]
        for f,m in fl:
            sg=sign_on_cube(f); info.append((f.degrees(),m,sg))
            if sg is None: s=None; break
            if m%2: s*=sg
        print(w,nm,'const',c,'factors',info,'=> sign',s,flush=True)

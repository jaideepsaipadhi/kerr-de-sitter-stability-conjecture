# Adaptive multi-chart driver: cover the (T, zeta, y) unit cube by cells; on each cell some chart must certify over its
# whole p-range.  Every chart parametrises a SUPERSET of the feasible set, so any one certifying chart suffices.
import pickle, sympy as sp, sys, random, time, json
from bern_core import *
T,ze,y,p=sp.symbols('T ze y p',positive=True); V=(T,ze,y,p)
def load(name):
    res=pickle.load(open(f'surf_{name}.pkl','rb')); sg={}
    for k,(P,dn) in res.items():
        vals=set(int(sp.sign(dn.subs({v:sp.Rational(random.randint(1,999),1000) for v in V}))) if hasattr(dn,'subs') else int(sp.sign(dn)) for _ in range(12))
        assert len(vals)==1,(name,k); sg[k]=vals.pop()
    BB={k:to_bern(sp.Poly(sp.expand(P.as_expr()*sg[k]),*V),V) for k,(P,dn) in res.items() if not P.is_zero}
    return BB
def rule(b):
    for k in ('bulk','c2','floor','Apos'):
        if k in b and mx(b[k])<0: return True
    if 'pp' in b and mx(b['pp'])<=0: return True
    return (mx(b['mup'])<=0 or mx(b['pos'])<=0) and (mn(b['mup'])>=0 or mx(b['neg'])<=0)
def restrict(BB,cell):
    out={}
    for k,A in BB.items():
        for d,(a,c) in enumerate(cell):
            # cut [a,c] out of [0,1] in dimension d by two splits at fractions (exact dyadic cells only)
            pass
        out[k]=A
    return out
def sub(BB,path):
    # path: list of (dim, side) dyadic subdivision steps
    out=BB
    for d,sd in path:
        out={k:split(A,d)[sd] for k,A in out.items()}
    return out
if __name__=='__main__':
    charts=[c for c in sys.argv[1].split(',')]; LOG='cert_'+charts[0]+'.log'; budget=float(sys.argv[2]); tmax=float(sys.argv[3])
    resume=len(sys.argv)>4 and sys.argv[4]=='resume'
    T0=time.time(); BBs={c:load(c) for c in charts}
    t0=time.time(); todo=json.load(open('left_'+charts[0]+'.json')) if resume else [[]]; done=0; fails=[]
    print('loaded',round(time.time()-T0,1),flush=True)
    while todo and time.time()-t0<tmax:
        path=todo.pop()
        okc=None
        for c in charts:
            r,deep=bb(sub(BBs[c],path),rule,4,tmax=budget,maxdep=40)
            if r['certified']: okc=c; break
        if okc:
            done+=1; open(LOG,'a').write(json.dumps([path,okc])+'\n'); continue
        if len(path)>=24: fails.append(path); continue
        # split the (T,ze,y) cell along the dimension cycled by depth
        d=len(path)%3
        todo.append(path+[(d,0)]); todo.append(path+[(d,1)])
    print(json.dumps({'complete':not todo and not fails,'cells_done':done,'todo':len(todo),'fails':len(fails),'sec':round(time.time()-t0,1)}))
    json.dump(todo+fails,open('left_'+charts[0]+'.json','w'))

import sympy as sp, pickle, sys, numpy as np, time
sys.path.insert(0,'/home/claude/kds/cert2')
from bern_core import to_bern, split, mn, mx, bb
r0,v=sp.symbols('r0 v',positive=True)
out=pickle.load(open('w2_poly.pkl','rb'))
for key,(c,facs) in out.items():
    small=[f for f,m in facs if sp.Poly(f,r0,v).total_degree()<=2]
    print(key,'small:',[ (sp.factor(f),m) for f,m in facs if sp.Poly(f,r0,v).total_degree()<=2])
    core=[f for f,m in facs if sp.Poly(f,r0,v).total_degree()>2][0]
    P=sp.Poly(core,r0,v); B=to_bern(P,[r0,v])
    print('   core Bernstein min/max sign:',np.sign(mn(B)),np.sign(mx(B)), 'corners:',[sp.sign(core.subs({r0:a,v:b})) for a in (0,1) for b in (0,1)], flush=True)

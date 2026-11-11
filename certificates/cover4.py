import flint, itertools, time, pickle, sympy as sp, sys
from fractions import Fraction as Fr
ctx=flint.fmpz_mpoly_ctx.get(('X','T','s','r'),'lex'); G=ctx.gens(); ONE=ctx.from_dict({(0,0,0,0):1})
def from_sympy(P):
    return ctx.from_dict({tuple(int(e) for e in m):int(c) for m,c in zip(P.monoms(),P.coeffs())})
def bern(p,deg):
    # coefficients (signs) of sum a_j x^j (1+x)^{n-j} per variable, with fixed degree box deg
    r=ctx.from_dict({tuple(n-int(e) for n,e in zip(deg,ex)):c for ex,c in p.to_dict().items()})
    s=r.compose(*[g+ONE for g in G],ctx=ctx)
    return {tuple(n-int(e) for n,e in zip(deg,ex)):int(c) for ex,c in s.to_dict().items()}
def half(p,k,right,deg):
    n=deg[k]; q=ctx.from_dict({ex:c*(2**(n-int(ex[k]))) for ex,c in p.to_dict().items()})   # 2^n p(x/2)
    if right: q=q.compose(*[(g+ONE) if i==k else g for i,g in enumerate(G)],ctx=ctx)       # 2^n p((x+1)/2)
    return q
def ok(p,deg,box):
    B=bern(p,deg)
    if min(B.values(),default=0)<0: return False
    for sel in itertools.product((0,1,2),repeat=4):
        if any((s==0 and box[k][0]==0) or (s==1 and box[k][1]==1) for k,s in enumerate(sel)): continue
        if not any(v>0 for ix,v in B.items() if all((s==2) or (s==0 and ix[k]==0) or (s==1 and ix[k]==deg[k]) for k,s in enumerate(sel))): return False
    return True
def cover(polys,degs,tmax=50,maxdep=40):
    # each box must satisfy ok() for at least one of polys
    t0=time.time(); st=[(polys,((Fr(0),Fr(1)),)*4,0)]; done=0
    while st:
        if time.time()-t0>tmax: return dict(ok=False,timeout=True,done=done,open=len(st))
        ps,box,dep=st.pop()
        if any(ok(p,d,box) for p,d in zip(ps,degs)): done+=1; continue
        if dep>=maxdep: return dict(ok=False,stuck=[(float(a),float(b)) for a,b in box])
        k=max(range(4),key=lambda i:box[i][1]-box[i][0]); a,c=box[k]; m=(a+c)/2
        for r,(lo,hi) in ((0,(a,m)),(1,(m,c))):
            nb=list(box); nb[k]=(lo,hi)
            st.append(([half(p,k,r,d) for p,d in zip(ps,degs)],tuple(nb),dep+1))
    return dict(ok=True,boxes=done,sec=round(time.time()-t0,1))
if __name__=='__main__':
    D=pickle.load(open('P173.pkl','rb'))
    P={k:from_sympy(v[0]) for k,v in D.items()}
    deg=lambda p: tuple(int(x) for x in p.degrees())
    which=sys.argv[1]
    if which=='g1': print('g1<0:',cover([-P['g1']],[deg(P['g1'])]))
    if which=='I':  print('I:',cover([P['g0'],P['q0']],[deg(P['g0']),deg(P['q0'])]))   # g0<0 <=> n_g0>0
    if which=='II': print('II:',cover([P['g0'],P['W']],[deg(P['g0']),deg(P['W'])]))
def at_r1(p):
    d={}
    for ex,c in p.to_dict().items():
        k=(int(ex[0]),int(ex[1]),int(ex[2]),0); d[k]=d.get(k,0)+int(c)
    return ctx.from_dict(d)
if __name__=='__main__' and len(sys.argv)>2 and sys.argv[2]=='r1':
    q0r=at_r1(P['q0']); Wr=at_r1(P['W'])
    print('I@rho=1:',cover([P['g0'],q0r],[deg(P['g0']),deg(q0r)]))
    print('II@rho=1:',cover([P['g0'],Wr],[deg(P['g0']),deg(Wr)]))

# Fast sign test: all Bernstein coeffs of p on [0,1]^3 >= 0 (some >0)  <=>  all coeffs of rev(shift(rev(p))) >= 0
import flint
def rev(p,ctx):
    d=p.degrees(); return ctx.from_dict({tuple(int(n)-int(e) for n,e in zip(d,ex)):c for ex,c in p.to_dict().items()})
def bern_signs(p):
    ctx=p.context(); g=ctx.gens(); one=ctx.from_dict({(0,)*len(g):1})
    r=rev(p,ctx); s=r.compose(*[x+one for x in g],ctx=ctx)
    # keep degree box: rev with original degrees
    d=p.degrees(); out=ctx.from_dict({tuple(int(n)-int(e) for n,e in zip(d,ex)):c for ex,c in s.to_dict().items()})
    cs=[int(c) for c in out.to_dict().values()]
    return min(cs), max(cs), len(cs)
def pos_open_cube(p):
    mn,mx,n=bern_signs(p); return mn>=0 and mx>0
def half(p,k,right):
    ctx=p.context(); g=list(ctx.gens()); one=ctx.from_dict({(0,)*len(g):1})
    if right: p=p.compose(*[(x+one) if i==k else x for i,x in enumerate(g)],ctx=ctx)
    n=int(p.degrees()[k])
    return ctx.from_dict({ex:c*(2**(n-int(ex[k]))) for ex,c in p.to_dict().items()})
def subdiv_pos(p,maxdepth=12,tmax=50):
    import time; t0=time.time(); st=[(p,0,(0,0,0),(0,0,0))]; boxes=0
    while st:
        if time.time()-t0>tmax: return None,boxes,'timeout'
        q,dep,lev,idx=st.pop()
        mn,mx,_=bern_signs(q)
        if mn>=0 and mx>0: boxes+=1; continue
        if mn<0 and _neg_vertex(q): return False,boxes,(lev,idx)
        if dep>=maxdepth: return None,boxes,('deep',lev,idx)
        k=min(range(3),key=lambda i:lev[i])
        for r in (0,1):
            l2=list(lev); l2[k]+=1; i2=list(idx); i2[k]=2*i2[k]+r
            st.append((half(q,k,r),dep+1,tuple(l2),tuple(i2)))
    return True,boxes,None
def _neg_vertex(q):
    # value at box corners (x in {0,1}^3): corner (0,0,0) is constant term; check all corners
    import itertools
    for c in itertools.product((0,1),repeat=3):
        v=sum(int(cf) for ex,cf in q.to_dict().items() if all(int(e)==0 or c[i]==1 for i,e in enumerate(ex)))
        if v<0: return True
    return False

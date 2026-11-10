# Exact certificate for p>0 on the open cube (x,v,tau) in (0,1)^3 for polynomials vanishing on cube edges.
# Near vanishing edges/corners: monomial domination (minimal Newton monomials with positive coefficient bounds
# dominate all other monomials on a small box); elsewhere: Bernstein B&B with de Casteljau.
import flint, numpy as np, time, sys, itertools
from fractions import Fraction
from math import comb, lcm
sys.path.insert(0,'/home/claude/kds/cert2'); from bern_core import split, mn, mx
c3=flint.fmpz_mpoly_ctx.get(('x','v','t'),'lex'); GX,GV,GT=c3.gens(); ONE=c3.from_dict({(0,0,0):1})
def local(p,flips):
    # flips[k]=1 means use s_k = 1 - coordinate_k ; return poly in (s0,s1,s2)
    sub=[(ONE-g) if f else g for g,f in zip((GX,GV,GT),flips)]
    return p.compose(*sub,ctx=c3)
def bern1(coeffs):  # coeffs list (monomial in u) -> exact Bernstein ints on [0,1] with common positive scale -> (min,max) as Fractions
    n=len(coeffs)-1
    B=[sum(Fraction(coeffs[j])*Fraction(comb(i,j),comb(n,j)) for j in range(i+1)) for i in range(n+1)]
    return min(B),max(B)
def poly_u(d,uidx,rng):
    # d: {exp: coeff} for monomials of a fixed transverse exponent; return coefficient list in u rescaled to rng
    lo,hi=rng; n=max([e[uidx] for e in d]+[0]); c=[Fraction(0)]*(n+1)
    for e,v in d.items(): c[e[uidx]]+=Fraction(int(v))
    # substitute u = lo + (hi-lo) w
    out=[Fraction(0)]*(n+1)
    for k,ck in enumerate(c):
        for j in range(k+1): out[j]+=ck*comb(k,j)*Fraction(lo)**(k-j)*Fraction(hi-lo)**j
    return out
def dominate(p,flips,trans,uidx,urng,eps):
    """trans: indices of transverse vars (s_i in [0,eps]); uidx: free var index or None (corner), urng its range."""
    q=local(p,flips).to_dict()
    groups={}
    for e,v in q.items():
        key=tuple(int(e[i]) for i in trans); groups.setdefault(key,{})[e]=v
    keys=list(groups); mins=[k for k in keys if not any(all(o[i]<=k[i] for i in range(len(k))) and o!=k for o in keys)]
    lb={}; ub={}
    for k,g in groups.items():
        if uidx is None:
            val=Fraction(int(sum(int(v) for v in g.values()))) ; lo_=hi_=val
        else:
            lo_,hi_=bern1(poly_u(g,uidx,urng))
        lb[k]=lo_; ub[k]=max(abs(lo_),abs(hi_))
    if any(lb[m]<=0 for m in mins): return False,'min coeff not positive'
    budget={m:lb[m] for m in mins}
    for k in keys:
        if k in mins: continue
        cands=[m for m in mins if all(m[i]<=k[i] for i in range(len(k)))]
        if not cands: return False,'undominated'
        best=max(cands,key=lambda m: sum(k[i]-m[i] for i in range(len(k))))
        extra=int(sum(int(k[i])-int(best[i]) for i in range(len(k))))
        budget[best]-=ub[k]*Fraction(eps)**extra
    ok=all(b>0 for b in budget.values())
    return ok, min(budget.values())
def dominate_adaptive(p,flips,trans,uidx,urng,eps,depth=0,maxd=14,stats=None):
    stats=stats if stats is not None else {'leaves':0}
    ok,_=dominate(p,flips,trans,uidx,urng,eps)
    if ok: stats['leaves']+=1; return True,stats
    if depth>=maxd: stats['fail']=urng; return False,stats
    lo,hi=urng; m=(lo+hi)/2
    a,_=dominate_adaptive(p,flips,trans,uidx,(lo,m),eps,depth+1,maxd,stats)
    if not a: return False,stats
    return dominate_adaptive(p,flips,trans,uidx,(m,hi),eps,depth+1,maxd,stats)
def dom_tubes(p,flips,trans,uidx,urng,elist,depth=0,maxd=12,out=None):
    """Return list of (urng,eps) covering urng, each certified by dominate, or None."""
    out=[] if out is None else out
    for e in elist:
        if dominate(p,flips,trans,uidx,urng,e)[0]: out.append((urng,e)); return out
    if depth>=maxd: return None
    lo,hi=urng; m=(lo+hi)/2
    if dom_tubes(p,flips,trans,uidx,(lo,m),elist,depth+1,maxd,out) is None: return None
    return dom_tubes(p,flips,trans,uidx,(m,hi),elist,depth+1,maxd,out)
def dom_tubes2(p,flips,trans,uidx,urng,elist,maxd=8):
    todo=[urng]; out=[]
    for e in elist:
        nxt=[]
        for r in todo:
            st=[(r,0)]
            while st:
                rr,d=st.pop()
                if dominate(p,flips,trans,uidx,rr,e)[0]: out.append((rr,e)); continue
                if d>=maxd: nxt.append(rr); continue
                lo,hi=rr; m=(lo+hi)/2; st+= [((lo,m),d+1),((m,hi),d+1)]
        todo=nxt
        if not todo: return out
    return None
def dominate2(p,flips,trans,uidx,urng,epsv,_q={}):
    """anisotropic: s_trans[i] in [0,epsv[i]]; greedy charging to the cheapest minimal monomial."""
    key0=(id(p),flips)
    if key0 not in _q: _q[key0]=local(p,flips).to_dict()
    q=_q[key0]; groups={}
    for e,v in q.items():
        groups.setdefault(tuple(int(e[i]) for i in trans),{})[e]=v
    keys=list(groups); mins=[k for k in keys if not any(all(o[i]<=k[i] for i in range(len(k))) and o!=k for o in keys)]
    lb={};ub={}
    for k,g in groups.items():
        lo_,hi_=bern1(poly_u(g,uidx,urng)); lb[k]=lo_; ub[k]=max(abs(lo_),abs(hi_))
    if any(lb[m]<=0 for m in mins): return False,'minneg'
    budget={m:lb[m] for m in mins}
    def w(k,m):
        r=Fraction(1)
        for i in range(len(k)): r*=Fraction(epsv[i])**(k[i]-m[i])
        return r
    for k in keys:
        if k in mins or ub[k]==0: continue
        c=[m for m in mins if all(m[i]<=k[i] for i in range(len(k)))]
        if not c: return False,'undom'
        best=min(c,key=lambda m: ub[k]*w(k,m)/lb[m]); budget[best]-=ub[k]*w(k,best)
    return all(b>0 for b in budget.values()), min(budget.values())

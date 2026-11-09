# Exact Bernstein branch-and-bound for B3 on the A = 0 surface over the COMPACT box
#   T = 1/t in [0,1], ze in [0,1], y in [0,1], bt = floor/B in [0,1].
# Certify: num_pos >= 0 where mu' >= 0 possible, num_neg >= 0 where mu' <= 0 possible, unless a constraint
# (R3 c2, R3 pp, bulk w>=10) is certainly violated on the box.  All arithmetic exact (python ints).
import sympy as sp, pickle, numpy as np, time, sys, itertools
from fractions import Fraction
T,ze,y,be=sp.symbols('T ze y be',positive=True); bt=sp.Symbol('bt',positive=True)
V=(T,ze,y,bt)
def polyT(expr):
    nu,dn=sp.fraction(sp.factor(sp.together(expr)))
    return sp.Poly(sp.expand(nu),*V), sp.factor(dn)
t0=0
t=1/T; B=(sp.Rational(11,2)*t-sp.Rational(5,4))/bt
x0=ze*(t-1); z0=1+x0; z=1+y*x0; s=(1-y)*x0; mu0=z0*x0*(t-z0)
cons={}
# each constraint as  P > 0 (or >= 0) with P a polynomial; we record the polynomial and the sign of its denominator
for nm,ex in (('c2',6*ze**2*(1-T)-T), ('pp',2*ze*((B-3)*T*bt)-bt*T), ('bulk',8*(B-3*z0)*s**2-10*mu0),
              ('mup',-3*z**2+2*(1+t)*z-t)):
    P,dn=polyT(sp.expand(ex)); cons[nm]=(P,dn); print(nm,'deg',P.degree(T),P.degree(ze),P.degree(y),P.degree(bt),'den',dn)

from math import comb
def to_bern(P):
    deg=[P.degree(v) for v in V]
    A=np.zeros([d+1 for d in deg],dtype=object); A[...]=Fraction(0)
    for mon,c in zip(P.monoms(),P.coeffs()): A[mon]=Fraction(int(sp.numer(c)),int(sp.denom(c)))
    for k,n in enumerate(deg):
        A=np.moveaxis(A,k,0); Bc=np.empty_like(A)
        for i in range(n+1):
            acc=0
            for j in range(i+1): acc=acc+A[j]*Fraction(comb(i,j),comb(n,j))
            Bc[i]=acc
        A=np.moveaxis(Bc,0,k)
    # common denominator -> ints
    from math import lcm
    L=1
    for v in A.flat: L=lcm(L,v.denominator)
    return np.vectorize(lambda v:int(v*L),otypes=[object])(A)
def split(A,k):
    A=np.moveaxis(A,k,0); n=A.shape[0]-1
    left=[A[0]*(2**n)]; right=[None]*(n+1); cur=[A[i]*(2**n)//1 for i in range(n+1)]  # scaled copies
    cur=[A[i] for i in range(n+1)]; L_=[cur[0]]; R_=[cur[n]]
    # de Casteljau with exact halving: keep everything multiplied by 2^n at the end
    levels=[cur]
    for r in range(1,n+1):
        prev=levels[-1]; nxt=[prev[i]+prev[i+1] for i in range(len(prev)-1)]  # = 2 * average
        levels.append(nxt)
    # level r entries are 2^r times true values; rescale to common 2^n
    Lc=[levels[r][0]*(2**(n-r)) for r in range(n+1)]
    Rc=[levels[n-r][r]*(2**(n-(n-r))) for r in range(n+1)]
    Lc=np.stack(Lc); Rc=np.stack(Rc)
    return np.moveaxis(Lc,0,k), np.moveaxis(Rc,0,k)
def mn(A): return min(A.flat)
def mx(A): return max(A.flat)
if __name__=='__main__':
    out=pickle.load(open('a0_scaled.pkl','rb'))
    polys={}
    for nm in ('pos','neg'):
        nu,_=out[nm]; P=sp.Poly(nu,T,ze,y,be); db=P.degree(be)
        polys[nm]=sp.Poly(sp.expand(nu.subs(be,1/bt)*bt**db),*V)
    for nm in ('c2','pp','bulk','mup'): polys[nm]=cons[nm][0]
    BB={nm:to_bern(P) for nm,P in polys.items()}
    print('bernstein built',round(time.time()-t0 if 't0' in dir() else 0,1),flush=True)
    deep=[]; start=time.time(); stack=[(BB,0,((0,1),(0,1),(0,1),(0,1)))]; closed=0; maxdepth=0
    while stack and time.time()-start<48:
        bb,dep,bx=stack.pop(0) if False else stack.pop(); maxdepth=max(maxdepth,dep)
        if dep>60: deep.append(bx); continue
        infeas=mx(bb['c2'])<0 or mx(bb['pp'])<=0 or mx(bb['bulk'])<0
        okp=mx(bb['mup'])<=0 or mn(bb['pos'])>=0
        okn=mn(bb['mup'])>=0 or mn(bb['neg'])>=0
        if infeas or (okp and okn): closed+=1; continue
        k=dep%4
        parts=[split(bb[nm],k) for nm in bb]
        a_,b_=bx[k]; m_=(a_+b_)/2
        bl=list(bx); bl[k]=(a_,m_); br=list(bx); br[k]=(m_,b_)
        stack.append(({nm:p[0] for nm,p in zip(bb,parts)},dep+1,tuple(bl)))
        stack.append(({nm:p[1] for nm,p in zip(bb,parts)},dep+1,tuple(br)))
    print('deep boxes',len(deep),[tuple((round(a,6),round(b,6)) for a,b in d) for d in deep[:5]]); print({'certified':not stack,'closed':closed,'open':len(stack),'maxdepth':maxdepth,'sec':round(time.time()-start,1)})

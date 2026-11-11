# Minimal vectorised interval arithmetic with outward rounding (np.nextafter after every op).
import numpy as np
INF=np.inf
def _dn(a): return np.nextafter(a,-INF)
def _up(a): return np.nextafter(a,INF)
class I:
    __slots__=('lo','hi')
    def __init__(s,lo,hi=None):
        lo=np.asarray(lo,float); hi=lo if hi is None else np.asarray(hi,float)
        s.lo,s.hi=lo,hi
    @staticmethod
    def c(v): return I(v,v)
    def __add__(s,o):
        o=o if isinstance(o,I) else I.c(o); return I(_dn(s.lo+o.lo),_up(s.hi+o.hi))
    __radd__=__add__
    def __neg__(s): return I(-s.hi,-s.lo)
    def __sub__(s,o):
        o=o if isinstance(o,I) else I.c(o); return I(_dn(s.lo-o.hi),_up(s.hi-o.lo))
    def __rsub__(s,o): return I.c(o)-s
    def __mul__(s,o):
        o=o if isinstance(o,I) else I.c(o)
        p=np.stack([s.lo*o.lo,s.lo*o.hi,s.hi*o.lo,s.hi*o.hi])
        bad=np.isnan(p).any(0)
        lo=np.where(bad,-INF,_dn(np.nanmin(np.where(np.isnan(p),INF,p),0)))
        hi=np.where(bad,INF,_up(np.nanmax(np.where(np.isnan(p),-INF,p),0)))
        return I(lo,hi)
    __rmul__=__mul__
    def inv(s):
        bad=(s.lo<=0)&(s.hi>=0)
        lo=np.where(bad,-INF,_dn(1/s.hi)); hi=np.where(bad,INF,_up(1/s.lo))
        return I(lo,hi)
    def __truediv__(s,o):
        o=o if isinstance(o,I) else I.c(o); return s*o.inv()
    def __rtruediv__(s,o): return I.c(o)*s.inv()
    def sq(s):
        a=np.abs(s.lo); b=np.abs(s.hi); m=np.maximum(a,b)
        lo=np.where((s.lo<=0)&(s.hi>=0),0.0,np.minimum(a,b)**2)
        return I(_dn(lo),_up(m*m))

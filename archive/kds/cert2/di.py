# Forward-mode AD over intervals: value interval + gradient intervals (for mean-value forms).
import numpy as np
from iv import I
def _I(c): return c if isinstance(c,I) else I.c(c)
class DI:
    __slots__=('v','g')
    def __init__(s,v,g): s.v=v; s.g=g
    @staticmethod
    def const(c,n):
        z=I.c(np.zeros_like(np.asarray(c,float))) if not isinstance(c,I) else I.c(np.zeros_like(c.lo))
        return DI(_I(c),[z]*n)
    def _w(s,o): return o if isinstance(o,DI) else DI.const(o,len(s.g))
    def __add__(s,o): o=s._w(o); return DI(s.v+o.v,[a+b for a,b in zip(s.g,o.g)])
    __radd__=__add__
    def __neg__(s): return DI(-s.v,[-a for a in s.g])
    def __sub__(s,o): o=s._w(o); return DI(s.v-o.v,[a-b for a,b in zip(s.g,o.g)])
    def __rsub__(s,o): return s._w(o)-s
    def __mul__(s,o): o=s._w(o); return DI(s.v*o.v,[a*o.v+s.v*b for a,b in zip(s.g,o.g)])
    __rmul__=__mul__
    def inv(s):
        iv=s.v.inv(); f=-(iv.sq()); return DI(iv,[f*a for a in s.g])
    def __truediv__(s,o): o=s._w(o); return s*o.inv()
    def __rtruediv__(s,o): return s._w(o)*s.inv()
    def sq(s): return DI(s.v.sq(),[2*s.v*a for a in s.g])
def hull(a,b):
    lo=np.minimum(a.lo,b.lo); hi=np.maximum(a.hi,b.hi); return I(lo,hi)
def dhull(a,b): return DI(hull(a.v,b.v),[hull(x,y) for x,y in zip(a.g,b.g)])

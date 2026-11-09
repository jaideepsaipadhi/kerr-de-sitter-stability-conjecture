import sympy as sp
class R:
    __slots__=('n','d')
    def __init__(s,n,d=None,V=None):
        s.n=n; s.d=d if d is not None else sp.Poly(1,*n.gens)
    @staticmethod
    def c(val,V): return R(sp.Poly(val,*V))
    def _w(s,o): return o if isinstance(o,R) else R(sp.Poly(o,*s.n.gens))
    def __add__(s,o):
        o=s._w(o)
        if s.d==o.d: return R(s.n+o.n,s.d)
        return R(s.n*o.d+o.n*s.d,s.d*o.d)
    __radd__=__add__
    def __neg__(s): return R(-s.n,s.d)
    def __sub__(s,o): return s+(-s._w(o))
    def __rsub__(s,o): return s._w(o)-s
    def __mul__(s,o): o=s._w(o); return R(s.n*o.n,s.d*o.d)
    __rmul__=__mul__
    def inv(s): return R(s.d,s.n)
    def __truediv__(s,o): o=s._w(o); return R(s.n*o.d,s.d*o.n)
    def __rtruediv__(s,o): return s._w(o)/s
    def __pow__(s,k): return R(s.n**k,s.d**k)

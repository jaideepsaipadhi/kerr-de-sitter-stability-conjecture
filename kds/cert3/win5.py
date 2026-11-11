import flint, time, pickle
t0=time.time()
ctx=flint.fmpz_mpoly_ctx.get(('r0','r2','b','tau','sg'),'lex')
r0,r2,b,tau,sg=ctx.gens()
one=ctx.from_dict({(0,0,0,0,0):1})
class F:
    __slots__=('n','d')
    def __init__(s,n,d=None):
        if isinstance(n,int): n=one*n
        if d is None: d=one
        g=n.gcd(d)
        if g!=one: n=n//g if (n%g==0) else n; d=d//g
        s.n=n; s.d=d
    def __add__(s,o):
        o=o if isinstance(o,F) else F(o)
        if s.d==o.d: return F(s.n+o.n,s.d)
        g=s.d.gcd(o.d); return F(s.n*(o.d//g)+o.n*(s.d//g), s.d*(o.d//g))
    __radd__=__add__
    def __neg__(s): return F(-s.n,s.d)
    def __sub__(s,o): return s+(-(o if isinstance(o,F) else F(o)))
    def __rsub__(s,o): return F(o)-s
    def __mul__(s,o):
        o=o if isinstance(o,F) else F(o)
        g1=s.n.gcd(o.d); g2=o.n.gcd(s.d)
        return F((s.n//g1)*(o.n//g2),(s.d//g2)*(o.d//g1))
    __rmul__=__mul__
    def inv(s): return F(s.d,s.n)
    def __truediv__(s,o): o=o if isinstance(o,F) else F(o); return s*o.inv()
    def __rtruediv__(s,o): return F(o)*s.inv()
r1=one; r3=-(r0+r1+r2); rs=[r0,r1,r2,r3]
e2=sum((rs[i]*rs[j] for i in range(4) for j in range(i+1,4)),0*one); e4=r0*r1*r2*r3
def prod(l):
    p=one
    for x in l: p=p*x
    return p
Pi=[prod([rs[j]-rs[k] for k in range(4) if k!=j]) for j in range(4)]
Pa=[Pi[0],-Pi[1],Pi[2],-Pi[3]]
rho=[x*x+b for x in rs]
y2=F(Pa[0]-Pa[1], rho[1]*Pa[0]-rho[0]*Pa[1]); y1=F(Pa[3]+Pa[1], rho[1]*Pa[3]+rho[3]*Pa[1])
y=y2+F(tau)*(y1-y2)
U=[(1-F(rho[0])*y)/F(Pa[0]),(F(rho[1])*y-1)/F(Pa[1]),(1-F(rho[2])*y)/F(Pa[2]),(F(rho[3])*y-1)/F(Pa[3])]
f1=U[1]+U[3]; f2=U[2]+U[3]
t=F((r1-r3)*(r2-r0),(r1-r0)*(r2-r3))
R=(f1+t*f2)/(f1+f2); z0=t/R
print('R z0',round(time.time()-t0,1),flush=True)
Xi=1-F(b*b,e4); al=F(-b*b,e4); lamv=F(-b,e4); K=F(one,(r1-r0)*(r2-r3))
LQ2=4*Xi*Xi*F(b)*U[0]*U[3]*R/lamv + K*(1-2*Xi*Xi*F(b)*y)
print('LQ2',round(time.time()-t0,1),flush=True)
LQ1=K*(1-4*al*Xi-4*Xi*Xi*F(b)*y)
LQ0=-K*(1+7*al)+lamv*(1+t-3*z0)
ex={'Q2':LQ2,'H2':LQ2+LQ1*F(1,2*one)+LQ0*F(1,4*one),'H1':LQ2+LQ1+LQ0+4*K}
out={}
for k,e in ex.items():
    out[k]=(str(e.n),str(e.d)); print(k,'deg n',e.n.degrees(),'terms',len(e.n),round(time.time()-t0,1),flush=True)
pickle.dump((out,str(e2),str(e4)),open('win_flint.pkl','wb'))

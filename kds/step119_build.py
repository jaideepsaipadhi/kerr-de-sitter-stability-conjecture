import pickle, sympy as sp, time
t_,z0_,A_,B_=sp.symbols('t z0 A B',positive=True); T,ze=sp.symbols('T ze',positive=True)
Tj=pickle.load(open('/home/claude/kds/T118.pkl','rb'))
t=1/T; z0=1+ze*(t-1); x0=z0-1
Bf=sp.Rational(11,2)*t-sp.Rational(5,4); Ac2=z0*((t-1)/(2*x0**2)-3); P=3+(t-1)/(2*x0)
def coefs(j):
    Pj=sp.Poly(Tj[j],A_,B_); return Pj.coeff_monomial(1),Pj.coeff_monomial(A_),Pj.coeff_monomial(B_)
out={}
def rec(name,ex):
    nu,dn=sp.fraction(sp.factor(sp.together(ex)))
    out[name]=(sp.Poly(sp.expand(nu),T,ze),dn); print(name,'deg',out[name][0].degree(T),out[name][0].degree(ze),'den',str(dn)[:80],flush=True)
sub=lambda e,a0: e.subs({t_:t,z0_:z0}).subs({A_:a0,B_:Bf})
rec('Ac2',Ac2)                       # case gate: Ac2>=0 (case II) or <=0 (case I)
rec('G0_I',0-(P-Bf)); rec('G0_II',Ac2-(P-Bf))
for case,a0 in (('I',0),('II',Ac2)):
    al={j:sub(Tj[j],a0) for j in (2,3,4,5,6)}
    for j in (4,5,6): rec(f'al{j}_{case}',al[j])
    _,b2,g2=coefs(2); _,b3,g3=coefs(3)
    b2,g2,b3,g3=[c.subs({t_:t,z0_:z0}) for c in (b2,g2,b3,g3)]
    rec(f'al3_{case}',al[3]); rec(f'al2_{case}',al[2])
    rec(f'v2_{case}',al[2]*(-b3)-(-b2)*al[3]); rec(f'v3_{case}',al[2]*(-g3)-(-g2)*al[3])
pickle.dump(out,open('/home/claude/kds/p1_2d.pkl','wb'))

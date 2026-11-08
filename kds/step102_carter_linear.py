# Item 102 / 98(c): Lambda-vacuum on Carter's ansatz with ARBITRARY Delta_r(R), Delta_x(I).
import sympy as sp, time
t0=time.time()
t,psi,R,I,L=sp.symbols('t psi R I Lambda')
Dr=sp.Function('Dr')(R); Dx=sp.Function('Dx')(I)
S=R**2+I**2
X=[t,psi,R,I]
g=sp.zeros(4)
# ds^2 = S(dR^2/Dr + dI^2/Dx) + (Dx (dt + R^2 dpsi)^2 - Dr (dt - I^2 dpsi)^2)/S
w1=[1,-I**2,0,0]; w2=[1,R**2,0,0]
for i in range(4):
    for j in range(4):
        g[i,j]=(Dx*w2[i]*w2[j]-Dr*w1[i]*w1[j])/S
g[2,2]=S/Dr; g[3,3]=S/Dx
gi=sp.simplify(g.inv())
Gam=[[[sp.cancel(sum(gi[k,q]*(sp.diff(g[q,i],X[j])+sp.diff(g[q,j],X[i])-sp.diff(g[i,j],X[q])) for q in range(4))/2) for j in range(4)] for i in range(4)] for k in range(4)]
def Ric(i,j):
    return sp.cancel(sum(sp.diff(Gam[k][i][j],X[k])-sp.diff(Gam[k][i][k],X[j])+sum(Gam[k][k][q]*Gam[q][i][j]-Gam[k][j][q]*Gam[q][i][k] for q in range(4)) for k in range(4)))

Rm=sp.zeros(4)
Rl={}
for i in range(4):
    for j in range(i,4):
        Rl[(i,j)]=Ric(i,j); Rl[(j,i)]=Rl[(i,j)]
c=sp.Symbol('c')
ok=True
for a_ in range(4):
    for b_ in range(4):
        mix=sp.cancel(sum(gi[a_,k]*Rl[(k,b_)] for k in range(4)))
        if mix==0: continue
        es=mix.subs({Dr:c*Dr,Dx:c*Dx}).doit()
        d2=sp.simplify(sp.diff(es,c,2))
        ok = ok and d2==0
print('mixed Ricci R^a_b affine (linear) in (Dr,Dx):', ok, round(time.time()-t0,1), flush=True)
er,m,n,gm=sp.symbols('eps m n gamma')
sub={Dr:-L/3*R**4+er*R**2-2*m*R+gm, Dx:-L/3*I**4-er*I**2+2*n*I+gm}
res=[sp.simplify((sum(gi[a_,k]*Rl[(k,b_)] for k in range(4))-(L if a_==b_ else 0)).subs(sub).doit()) for a_ in range(4) for b_ in range(4)]
print('quartic family solves R^a_b = Lambda delta:', all(r==0 for r in res), round(time.time()-t0,1))

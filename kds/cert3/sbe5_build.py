import flint,time,pickle
exec(open("smu5.py").read().split("\nex={")[0])
exec("def bcoeffs"+open('margin_red5.py').read().split("def bcoeffs")[1].split("res={}")[0])
Sg=F(sg); T=t-z0; aa=z0-1; two=F(2*one)
gam=F(4*one)*U[0]*U[3]*z0*z0/(t*(f1+f2)*(f1+f2))     # = Gamma/T
alp=aa/T; io=F(one)/T
den=F(4*one)*Xi*Xi*F(b)*f2*(f1+t*f2)
Hm={1:LQ2p+LQ1+LQ0p+4*K, 2:LQ2p+LQ1*F(1,2*one)+LQ0p*F(1,4*one)}
P=two*alp*(alp+Sg)+(two*alp+Sg)*io
E=P*(1-Sg)-Sg*(alp+Sg)*(alp+Sg+io)
ex={}
for m,H in Hm.items():
    rho=F(21,10*one)*(-H)*lamv/den
    W1=gam*(alp+Sg)*(alp+Sg)-two*rho*Sg*P
    ex[f'R{m}']=-H; ex[f'W1_{m}']=W1; ex[f'W2_{m}']=W1-rho*rho*Sg*E
    ex[f'Z_{m}']=gam*T-rho*(F(27,10*one)+F(54,10*one)*aa)
res={}
for k,e in ex.items():
    cd,fld=e.d.factor(); cn,fln=e.n.factor(); sgn0=(1 if cn>0 else -1)*(1 if cd>0 else -1)
    facs=[]
    for f,mm in list(fln)+list(fld):
        if mm%2==0: continue
        if f.degrees()[2]==0: facs.append(('plain',f.to_dict()))
        else:
            F1c,F0c=reduce_b(f); Y=-(2*F0c+F1c*e2); D=e2*e2-4*e4; Pp=F1c*F1c*D-Y*Y
            facs.append(('b',(F1c.to_dict(),Y.to_dict(),Pp.to_dict(),None,None,None,f.to_dict())))
    res[k]=(sgn0,facs); print(k,[(fa[0]) for fa in facs],max(len(ctx.from_dict(fa[1][2] if fa[0]=='b' else fa[1])) for fa in facs),round(time.time()-t0,1),flush=True)
pickle.dump(res,open('margin/sbe5_res.pkl','wb'))

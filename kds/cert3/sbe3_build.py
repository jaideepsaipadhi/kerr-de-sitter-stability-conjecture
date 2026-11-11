import flint,time,pickle
exec(open("smu5.py").read().split("\nex={")[0])
exec("def bcoeffs"+open('margin_red5.py').read().split("def bcoeffs")[1].split("res={}")[0])
Sg=F(sg); T=t-z0; aa=z0-1; two=F(2*one)
Gam=F(4*one)*U[0]*U[3]*T*z0*z0/(t*(f1+f2)*(f1+f2))
den=F(4*one)*Xi*Xi*F(b)*f2*(f1+t*f2)
Hm={1:LQ2p+LQ1+LQ0p+4*K, 2:LQ2p+LQ1*F(1,2*one)+LQ0p*F(1,4*one)}
s_=T*Sg; z=z0+s_; c=aa+s_
Ep=(two*aa*z+s_)*(1-Sg)-Sg*z*c
ex={}
for m,H in Hm.items():
    rho=F(21,10*one)*(-H)*lamv/den
    W1=Gam*c*c-two*rho*s_*(two*aa*z+s_)
    W2=W1-F(7,10*one)*rho*s_*Ep
    ex[f'R{m}']=-H          # rho_m<=0  <=>  -H<=0 (den>0: f2<0,f1+tf2<0 -> f2(f1+tf2)>0)
    ex[f'W1_{m}']=W1; ex[f'W2_{m}']=W2
res={}
for k,e in ex.items():
    cd,fld=e.d.factor(); cn,fln=e.n.factor(); sgn0=(1 if cn>0 else -1)*(1 if cd>0 else -1)
    facs=[]
    for f,mm in list(fln)+list(fld):
        if mm%2==0: continue
        if f.degrees()[2]==0: facs.append(('plain',f.to_dict()))
        else:
            F1c,F0c=reduce_b(f); Y=-(2*F0c+F1c*e2); D=e2*e2-4*e4; P=F1c*F1c*D-Y*Y
            F1a=-F0c; Ya=-(2*e4*F1c+F0c*e2); Pa=F1a*F1a*D-Ya*Ya
            facs.append(('b',(F1c.to_dict(),Y.to_dict(),P.to_dict(),F1a.to_dict(),Ya.to_dict(),Pa.to_dict())))
    res[k]=(sgn0,facs); print(k,[f[0] for f in facs],round(time.time()-t0,1),flush=True)
pickle.dump(res,open('margin/sbe3_res.pkl','wb'))

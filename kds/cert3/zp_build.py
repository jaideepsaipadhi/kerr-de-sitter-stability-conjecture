exec(open("sal_build.py").read().split("\nex={")[0])
exec("def bcoeffs"+open('margin_red.py').read().split("def bcoeffs")[1].split("res={}")[0])
T=t-z0; aa=z0-1
gam=F(4*one)*U[0]*U[3]*z0*z0/(t*(f1+f2)*(f1+f2)); alp=aa/T; io=F(one)/T
den=F(4*one)*Xi*Xi*F(b)*f2*(f1+t*f2)
Hm={1:LQ2p+LQ1+LQ0p+4*K, 2:LQ2p+LQ1*F(1,2*one)+LQ0p*F(1,4*one)}
ex={}
for m,H in Hm.items():
    rho=F(21,10*one)*(-H)*lamv/den
    ex[f'Zp{m}']=gam*(1+alp)-(F(4*one)*rho+F(2*one)*rho*rho)*(alp+io); ex[f'NR{m}']=H*io
res={}
for k,e in ex.items():
    cd,fld=e.d.factor(); cn,fln=e.n.factor(); sgn0=(1 if cn>0 else -1)*(1 if cd>0 else -1)
    facs=[]
    for f,m in list(fln)+list(fld):
        if m%2==0: continue
        if f.degrees()[2]==0: facs.append(('plain',f.to_dict()))
        else:
            F1c,F0c=reduce_b(f); Y=-(2*F0c+F1c*e2); D=e2*e2-4*e4; P=F1c*F1c*D-Y*Y
            facs.append(('b',(F1c.to_dict(),Y.to_dict(),P.to_dict())))
        print(k,f.degrees(),m,facs[-1][0],flush=True)
    res[k]=(sgn0,facs)
pickle.dump(res,open('margin/zp_res.pkl','wb'))

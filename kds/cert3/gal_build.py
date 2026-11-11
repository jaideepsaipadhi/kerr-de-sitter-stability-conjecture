exec(open("sal_build.py").read().split("\nex={")[0])
exec("def bcoeffs"+open('margin_red.py').read().split("def bcoeffs")[1].split("res={}")[0])
# F per m per a: Fh ; F^2 = b*Fh^2*m^2 ; lamv*G_alpha/m^2-structure:
Bq2=K*(1-2*Xi*Xi*F(b)*y)
P=F(3*one)*(t-1)*F(b)*Fh*Fh*lamv
c=F(84,10*one)*t*(F(4*one)*t-2)
ex={'Q2':P+c*Bq2,'H2':P+c*(Bq2+LQ1*F(1,2*one)+LQ0p*F(1,4*one)),'H1':P+c*(Bq2+LQ1+LQ0p+4*K)}
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
pickle.dump(res,open('margin/gal_res.pkl','wb'))

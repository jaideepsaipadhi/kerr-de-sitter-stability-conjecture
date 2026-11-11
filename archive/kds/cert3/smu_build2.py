exec(open("smu_build.py").read().split("\nexec(\"def bcoeffs")[0])
exec("def bcoeffs"+open('margin_red.py').read().split("def bcoeffs")[1].split("res={}")[0])
ex={'LH2':ex['H2'],'LH1':ex['H1']}
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
pickle.dump(res,open('margin/smu2_res.pkl','wb'))

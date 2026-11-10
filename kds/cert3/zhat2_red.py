import flint, time, pickle, sys
exec(open('zhat2_build.py').read())
src=open('win_red3.py').read()
exec("def bcoeffs"+src.split("def bcoeffs")[1].split("res={}")[0])
t0=time.time()
for k in sys.argv[1:]:
    e={'q2':q2,'Hhalf':Hhalf,'H1':H1,'q0':q0}[k]
    F1c,F0c=reduce_b(e.n*e.d)
    Y=-(2*F0c+F1c*e2); D=e2*e2-4*e4; Pp=F1c*F1c*D-Y*Y
    print(k,'F1',F1c.degrees(),'P',Pp.degrees(),len(Pp),round(time.time()-t0,1),flush=True)
    pickle.dump({k:tuple(x.to_dict() for x in (F1c,Y,Pp))},open('zhat2_res_%s.pkl'%k,'wb'))

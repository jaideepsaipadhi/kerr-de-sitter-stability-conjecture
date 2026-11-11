import flint, time, pickle, sys
exec(open('zhat_build.py').read().split("out={k:")[0])
exec(open('win_red3.py').read().split("def bcoeffs")[1].join(["def bcoeffs",""]) if False else "")
src=open('win_red3.py').read()
exec("def bcoeffs"+src.split("def bcoeffs")[1].split("res={}")[0])
t0=time.time(); res={}
for k in sys.argv[1:]:
    e={'q2':q2,'Hhalf':Hhalf,'H1':H1}[k]
    F1c,F0c=reduce_b(e.n*e.d)
    Y=-(2*F0c+F1c*e2); D=e2*e2-4*e4; Pp=F1c*F1c*D-Y*Y
    res[k]=(F1c,Y,Pp); print(k,'F1',F1c.degrees(),len(F1c),'P',Pp.degrees(),len(Pp),round(time.time()-t0,1),flush=True)
    pickle.dump({kk:tuple(x.to_dict() for x in v) for kk,v in res.items()},open('zhat_res_%s.pkl'%k,'wb'))

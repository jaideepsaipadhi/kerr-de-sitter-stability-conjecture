import flint, time, pickle, sys
ctx=flint.fmpz_mpoly_ctx.get(('r0','r2','b','tau'),'lex'); r0,r2,b,tau=ctx.gens(); one=ctx.from_dict({(0,0,0,0):1})
r1=one; r3=-(r0+r1+r2); rs=[r0,r1,r2,r3]
e2=sum((rs[i]*rs[j] for i in range(4) for j in range(i+1,4)),0*one); e4=r0*r1*r2*r3
src=open('win_red3.py').read()
exec("def bcoeffs"+src.split("def bcoeffs")[1].split("res={}")[0])
E=pickle.load(open('zhat2_exprs.pkl','rb')); t0=time.time()
k=sys.argv[1]; n,d=[ctx.from_dict(x) for x in E[k]]
F1c,F0c=reduce_b(n*d); print('reduced',round(time.time()-t0,1),flush=True)
Y=-(2*F0c+F1c*e2)
pickle.dump({k:(F1c.to_dict(),Y.to_dict(),{})},open('zhat2_res_%s.pkl'%k,'wb')); print('saved',F1c.degrees(),round(time.time()-t0,1))

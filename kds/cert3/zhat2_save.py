import pickle
exec(open('zhat2_build.py').read())
pickle.dump({k:(v.n.to_dict(),v.d.to_dict()) for k,v in dict(q2=q2,q0=q0,Hhalf=Hhalf,H1=H1).items()},open('zhat2_exprs.pkl','wb'))
print('saved')

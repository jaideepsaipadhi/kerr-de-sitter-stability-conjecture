exec(open('step141b.py').read().split("for c in (0.5")[0])
fams={}
for p in (1.5,1.8,2.0):
  for c in (1.5,2.0,2.5,3.0):
    fams[('poly',p,c)]=(lambda x,p=p,c=c: x**p*(1-c*x*(1-x))*(1-c*x*(1-x)>0)+1e-300, lambda x,p=p,c=c: p*x**(p-1)*(1-c*x*(1-x))+x**p*(-c*(1-2*x)))
for a in (0.5,1,2):
    fams[('exp',2,a)]=(lambda x,a=a: x**2*np.exp(-a*x*(1-x)), lambda x,a=a: np.exp(-a*x*(1-x))*(2*x-a*x**2*(1-2*x)))
res={}
for k,(f,df) in fams.items():
    res[k]=[cond(Q,f,df) for Q in out]
anyok=np.any(np.array(list(res.values())),axis=0).sum()
for k,v in sorted(res.items(),key=lambda kv:-sum(kv[1]))[:6]: print(k,sum(v))
print('any',anyok,len(out))

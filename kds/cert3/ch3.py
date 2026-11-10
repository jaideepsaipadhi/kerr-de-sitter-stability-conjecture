from ch2 import *
F1,Y,P=res['Q2']; p=core3(F1)[0]; q,k=chart(p,'B'); f=fe(q)
g=np.linspace(0,1,21); best=[]
for a in g:
  for b_ in g:
    for c in g:
      val=f(a,b_,c)
      if val<0: best.append((a,b_,c,val))
print(len(best)); print(best[:10])
fp=fe(p)
for (a,b_,c,_) in best[:5]:
    if c>0: print('orig point x,v,tau',a,c,1-b_*c,'orig val',fp(a,c,1-b_*c))

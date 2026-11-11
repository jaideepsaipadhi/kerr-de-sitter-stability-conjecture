from ch2 import *
c1=flint.fmpz_mpoly_ctx.get(('s',),'lex'); (s,)=c1.gens()
for which in ('Q2','H2','H1'):
  F1,Y,P=res[which]
  for nm,d in (('F1',F1),('P',P)):
    for p in core3(d):
        one=c1.from_dict({(0,):1})
        e1=p.compose(one*0+one, s, one, ctx=c1)       # x=1,tau=1
        e2=p.compose(s, one*0, one, ctx=c1)            # v=0,tau=1
        e3=p.compose(one, one, s, ctx=c1)              # x=1,v=1 (tau free)
        e4=p.compose(s, one, one*0, ctx=c1)            # v=1,tau=0
        e5=p.compose(one, s, one*0, ctx=c1)            # x=1,tau=0
        print(which,nm,p.degrees(),'x=1,t=1:',e1==0,' v=0,t=1:',e2==0,' x=1,v=1:',e3==0,' v=1,t=0:',e4==0,' x=1,t=0:',e5==0)

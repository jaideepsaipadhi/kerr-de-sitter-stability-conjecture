from ch2 import *
import sys
which,nm,idx=sys.argv[1],sys.argv[2],int(sys.argv[3])
F1,Y,P=res[which]; p=core3(F1 if nm=='F1' else P)[idx]
c=flint.fmpz_mpoly_ctx.get(('X','v','T'),'lex'); X_,v_,T_=c.gens(); one=c.from_dict({(0,0,0):1})
q=p.compose(one-X_,v_,one-T_,ctx=c); d=q.to_dict(); pts=list(d)
low=[pt for pt in pts if not any(all(o[i]<=pt[i] for i in range(3)) and o!=pt for o in pts)]
print('minimal (X,v,T) exponents with coeff:',sorted((pt,int(d[pt])) for pt in low))

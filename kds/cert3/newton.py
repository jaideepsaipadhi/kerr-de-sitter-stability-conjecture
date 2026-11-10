from ch2 import *
import sys
which,nm,idx=sys.argv[1],sys.argv[2],int(sys.argv[3])
F1,Y,P=res[which]; p=core3(F1 if nm=='F1' else P)[idx]
cT=flint.fmpz_mpoly_ctx.get(('x','v','T'),'lex'); x_,v_,T_=cT.gens(); one=cT.from_dict({(0,0,0):1})
q=p.compose(x_,v_,one-T_,ctx=cT)   # tau = 1 - T
d=q.to_dict(); pts=sorted({(e[1],e[2]) for e in d})
# lower convex hull of exponents (j=v-power,k=T-power)
import itertools
low=[pt for pt in pts if not any(o[0]<=pt[0] and o[1]<=pt[1] and o!=pt for o in pts)]
print('minimal (v,T) exponents:',low)
for pt in low:
    coef=cT.from_dict({(e[0],0,0):c for e,c in d.items() if (e[1],e[2])==pt})
    vals=[float(coef.compose(cT.from_dict({(0,0,0):flint.fmpz(n)})*0+one*0+one*n if False else one*0, one*0, one*0, ctx=cT).to_dict().get((0,0,0),0)) for n in [0]]
    print(pt,'coef in x:',coef)

import sys,json; sys.argv=['x','margin/sal_res.pkl','Q2,H2,H1','sal','0','40']
exec(open('cover3d_f.py').read().split("PK=f'margin/c3f_")[0])
R=json.load(open('margin/c3f_sal.json')); D=[z for z in R if sum(z[0])>=36]
import random; random.seed(1)
for l,i in random.sample(D,4):
    ps=descend(tuple(l),tuple(i)); d=sum(l)
    print([round(i[k]/2**l[k],4) for k in range(3)],{k:[bsg(ps,(k,n),d) for n in ('F1','Y','P')] for k in TG}, [len(ps[(k,'Y')]) for k in TG])

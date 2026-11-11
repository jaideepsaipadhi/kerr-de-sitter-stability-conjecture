import sys,json; sys.argv=['x','margin/sal2_res.pkl','Q2,H2,H1','sal2','0','40']
exec(open('cover3d_g.py').read().split("PK=f'margin/c3f_")[0])
R=json.load(open('margin/c3f_sal2.json')); D=[z for z in R if sum(z[0])>=36]
import random; random.seed(2)
for l,i in random.sample(D,3):
    ps=descend(tuple(l),tuple(i)); d=sum(l)
    for k in ('Q2',):
        s0,its=ITEMS[k]
        print([round(i[q]/2**l[q],4) for q in range(3)],k,[(kind,j,bsg(ps,(k,j,'v'),d) if kind=='plain' else (bval(ps,k,j,d),[bsg(ps,(k,j,n),d) for n in ('F1','Y','P')])) for kind,j in its])

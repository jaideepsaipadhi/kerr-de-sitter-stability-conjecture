import sys,os,json; os.environ['ALT']='margin/smu2_res.pkl:LH2,LH1'; sys.argv=['x','margin/sal2_res.pkl','H2,H1','sal2','0','40','EB']
exec(open('cover3d_h.py').read().split("PK=f'margin/c3f_")[0])
R=json.load(open('margin/c3f_sal2EB.json')); D=[z for z in R if sum(z[0])>=27]
import random; random.seed(3)
for l,i in random.sample(D,3):
    ps=descend(tuple(l),tuple(i)); d=sum(l)
    for k in ('H2','H1'):
        s0,its=ITEMS[k]
        print([round(i[q]/2**l[q],4) for q in range(3)],k,[(j, bsg(ps,(k,j,'v'),d) if kind=='plain' else (bval(ps,(k,),j,d),[bsg(ps,(k,j,n),d) for n in ('F1','Y','P')])) for kind,j in its])

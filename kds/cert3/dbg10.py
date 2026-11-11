import sys,json; sys.argv=['x','margin/sbe2_res.pkl','sbe2','0','48']
exec(open('cover4d.py').read().split("PK=f'margin/c4_")[0])
R=json.load(open('margin/c4_sbe2.json')); R.sort(key=lambda z:-sum(z[0]))
for l,i in R[:3]:
    ps=descend(tuple(l),tuple(i)); d=sum(l)
    out={}
    for k in ITEMS:
        s0,its=ITEMS[k]
        out[k]=[(kind,bsg(ps,(k,j,'v'),d) if kind=='plain' else (bval(ps,k,j,d),[bsg(ps,(k,j,n),d) for n in ('F1','Y','P')])) for kind,j in its]
    print([round(i[q]/2**l[q],3) for q in range(4)],out)

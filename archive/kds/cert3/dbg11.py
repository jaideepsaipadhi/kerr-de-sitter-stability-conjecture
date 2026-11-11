import sys,json; sys.argv=['x','margin/sbe3_res.pkl','sbe3','0','48']
exec(open('cover4d.py').read().split("PK=f'margin/c4_")[0])
R=json.load(open('margin/c4_sbe3.json')); R.sort(key=lambda z:-sum(z[0]))
for l,i in R[:3]:
    ps=descend(tuple(l),tuple(i)); d=sum(l); BOX[0]=[(i[q]/2**l[q],(i[q]+1)/2**l[q]) for q in (0,1)]
    out={}
    for k in ITEMS:
        s0,its=ITEMS[k]
        out[k]=[(kind,bsg(ps,(k,j,'v'),d) if kind=='plain' else bval(ps,k,j,d)) for kind,j in its]
    print([round(i[q]/2**l[q],3) for q in range(4)],out)

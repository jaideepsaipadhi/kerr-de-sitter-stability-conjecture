import sys,os,json; os.environ['PERM2']='1'; os.environ['ALT']='margin/zp_res.pkl:NR1,NR2'
sys.argv=['x','margin/zp_res.pkl','Zp1,Zp2','zp','0','40','T1a']
exec(open('cover3d_h.py').read().split("PK=f'margin/c3f_")[0])
R=json.load(open('margin/c3f_zpT1a.json')); R.sort(key=lambda z:-sum(z[0]))
for l,i in R[:2]:
    ps=descend(tuple(l),tuple(i)); d=sum(l)
    for k in ('Zp2',):
        s0,its=ITEMS[k]
        print([round(i[q]/2**l[q],4) for q in range(3)],k,s0,[(kind,j,bsg(ps,(k,j,'v'),d) if kind=='plain' else (bval(ps,(k,),j,d),[bsg(ps,(k,j,n),d) for n in ('F1','Y','P')])) for kind,j in its])
    s0,its=ITEMS2['NR2'] if 'ITEMS2' in dir() else (None,[])

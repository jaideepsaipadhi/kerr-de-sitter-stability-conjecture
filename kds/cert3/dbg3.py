import sys,os; sys.argv=['x','EA','0','30']; os.environ['CONDS']='q0,q1,c2neg,v1,v2,disc,mm,pm,mp'
exec(open('cover2d_f.py').read().split("PK=f'margin/c2f_")[0])
for (L,I) in [((16,16),(int(0.176575*2**16),int(0.827545*2**16))),((16,16),(int(0.29920*2**16),int(0.05645*2**16)))]:
    ps=descend(L,I)
    for tk in ('Q2','H2','H1'):
        print(tk,{c:([bsg(ps,(tk,c,n),20) for n in ('F1','Y','P')],pos(ps,tk,c,20)) for c in ('q0','q1','c2neg','v1','v2','disc','mm','pm','mp')})

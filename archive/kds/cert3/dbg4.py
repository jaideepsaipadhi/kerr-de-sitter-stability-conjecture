import sys,os,json; sys.argv=['x','EA','0','30']; os.environ['CONDS']='q0,q1,c2neg,v1,v2,disc,mm,pm,mp'
exec(open('cover2d_f.py').read().split("PK=f'margin/c2f_")[0])
R=json.load(open('margin/c2f_EA.json')); l,i=max(R,key=lambda z:sum(z[0]))
for dl in range(0,8,2):
    L=(l[0]+dl,l[1]+dl); I=(i[0]<<dl,i[1]<<dl)
    ps=descend(L,I)
    for TK in ('H2','H1'): print(TK,L,{c:([bsg(ps,(TK,c,n),sum(L)) for n in ('F1','Y','P')],pos(ps,TK,c,sum(L))) for c in ('q0','q1','v1','c2neg','v2','disc','mm','pm','mp')})

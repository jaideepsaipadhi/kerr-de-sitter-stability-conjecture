import sys,os,json; sys.argv=['x','EA','0','30']; os.environ['CONDS']='q0,q1,c2neg,v1,v2,disc,mm,pm,mp'
exec(open('cover2d_f.py').read().split("PK=f'margin/c2f_")[0])
R=json.load(open('margin/c2f_EA.json'))
deep=[z for z in R if sum(z[0])>=24][:4]
for l,i in deep:
    ps=descend(tuple(l),tuple(i)); d=sum(l)
    for TK in ('Q2','H2','H1'):
        bad={c:[bsg(ps,(TK,c,n),d) for n in ('F1','Y','P')] for c in ('q0','q1','v1','c2neg','v2','disc','mm','pm','mp') if not pos(ps,TK,c,d)}
        print(l,[round(i[k]/2**l[k],6) for k in (0,1)],TK,bad)

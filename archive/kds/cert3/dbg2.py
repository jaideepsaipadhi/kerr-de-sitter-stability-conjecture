import sys,os; sys.argv=['x','EA','0','30']; os.environ['CONDS']='q0,q1,c2neg,v1,v2,disc'
src=open('cover2d_f.py').read().split("PK=f'margin/c2f_")[0]
exec(src)
ps=descend((6,6),(20,40))
for tk in ('H2',):
    for c in ('q0','q1','v1','disc','c2neg','v2'):
        print(c,[bsg(ps,(tk,c,n),12) for n in ('F1','Y','P')], [ [ (round(A.min(),3),round(A.max(),3)) for A in ps[(tk,c,n)]] for n in ('F1',)])
print('closed?',closed(ps,12))
for tk in ('Q2','H1'):
    for c in ('q0','q1','v1','disc'):
        print(tk,c,pos(ps,tk,c,12))

import sys,os,json; os.environ['PERM2']='1'; os.environ['ALT']='margin/zp_res.pkl:NR1,NR2'
sys.argv=['x','margin/zp_res.pkl','Zp1,Zp2','zp','0','40','T1a']
exec(open('cover3d_h.py').read().split("PK=f'margin/c3f_")[0])
R=json.load(open('margin/c3f_zpT1a.json')); R.sort(key=lambda z:-sum(z[0]))
l,i=R[0]; ps=descend(tuple(l),tuple(i)); d=sum(l)
BOX[0]=[(i[q]/2**l[q],(i[q]+1)/2**l[q]) for q in range(3)]
print('box',BOX[0])
print([k for k in S if k[0]=='Zp2' and k[1]==0])
print('bval0',bval0(ps,('Zp2',),0,d),'unred',unred(ps,('Zp2',),0,d))
bx=BOX[0]; print([xv_of((c0,c1,c2)) for c0 in bx[0] for c1 in bx[1] for c2 in bx[2]][:4])

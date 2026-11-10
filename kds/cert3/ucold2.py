import sys,pickle,numpy as np
sys.argv=['x','margin/sca_res.pkl','CQ2','scC','0','40','EA2']
exec(open('gcover.py').read().split("S={}")[0])
CH='EA2'
exec("def chart2"+open('gcover.py').read().split("def chart2")[1].split("root=")[0])
res=pickle.load(open('margin/sca_res.pkl','rb'))
def full(d):
    q=tocube(ctx4.from_dict(d)); 
    # chart without stripping s-power
    src=chart(q); return src
for k in ('CQ2',):
    F1=chart(tocube(ctx4.from_dict(res[k][0]))); Y=chart(tocube(ctx4.from_dict(res[k][1])))
    # chart() strips min s-power; recompute orders by comparing
    print('degrees',F1.degrees(),Y.degrees())

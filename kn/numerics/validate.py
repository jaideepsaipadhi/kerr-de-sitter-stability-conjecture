import sys; sys.argv=['x']
exec(open('realmode_scan.py').read().split("if __name__")[0])
av,Qv,mv=0.001,0.001,2; mats=setup(av,Qv,mv)
for W in [0.37367168-0.08896232j, 0.40-0.08896j, 0.37367-0.05j, 0.30-0.09j]:
    for rb in [100.,200.]:
        print(W, rb, f"{smin(av,Qv,mv,W,mats,rbig=rb,rmid=3.0):.3e}",flush=True)

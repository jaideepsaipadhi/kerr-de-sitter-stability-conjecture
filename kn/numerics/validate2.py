import sys; sys.argv=['x']
exec(open('realmode_scan.py').read().split("if __name__")[0])
av,Qv,mv=0.99,0.0001,2; mats=setup(av,Qv,mv)
for W in [0.8710-0.0294j,0.8710-0.02j,0.86-0.0294j,0.88-0.0294j,0.8710-0.04j]:
    print(W, f"{smin(av,Qv,mv,W,mats,rbig=40.,rmid=3.0,x0=1e-5):.3e}",flush=True)

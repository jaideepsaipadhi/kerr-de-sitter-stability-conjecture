exec(open("scg_build.py").read().split("\nres={}")[0])
for k,e in ex.items():
    print(k,'n deg',e.n.degrees(),'d deg',e.d.degrees())
    c,fl=e.d.factor(); print("const",c,[(f.degrees(),m) for f,m in fl])
    break

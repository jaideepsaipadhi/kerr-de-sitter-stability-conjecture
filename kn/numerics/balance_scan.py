import numpy as np, sys
exec(open('balance_kerr.py').read().split('# ---- near zone')[0].replace('s,m=-2,2','s,m=-2,int(sys.argv[1])'))
exec('# ---- near zone'+open('balance_kerr.py').read().split('# ---- near zone')[1].split('for n in')[0])
ns=[-50,-5,-1,-.1,-.01,-.001,.001,.1,1,5]
print(f"m={m}: |F|={abs(F):.4g}  |N|(n) =", "  ".join(f"{n}:{abs(near_ratio(n)):.4g}" for n in ns))

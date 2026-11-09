# verify that the logged dyadic cells (paths) tile the unit cube in (T, zeta, y)
import json, sys
from fractions import Fraction
cells=[json.loads(l)[0] for l in open(sys.argv[1])]
vol=Fraction(0)
boxes=[]
for path in cells:
    lo=[Fraction(0)]*3; hi=[Fraction(1)]*3
    for d,sd in path:
        m=(lo[d]+hi[d])/2
        if sd==0: hi[d]=m
        else: lo[d]=m
    boxes.append((lo,hi)); vol+=(hi[0]-lo[0])*(hi[1]-lo[1])*(hi[2]-lo[2])
# pairwise disjointness of interiors
def disj(a,b): return any(a[1][k]<=b[0][k] or b[1][k]<=a[0][k] for k in range(3))
ok=all(disj(boxes[i],boxes[j]) for i in range(len(boxes)) for j in range(i+1,len(boxes)))
print('cells',len(boxes),'total volume',vol,'disjoint',ok,'TILES CUBE' if vol==1 and ok else 'NOT A TILING')

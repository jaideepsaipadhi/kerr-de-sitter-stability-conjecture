import json,sys
R=json.load(open(sys.argv[1]))
def ge(z,k,a,b): return z[1][k]*b>=a*2**z[0][k]
def le(z,k,a,b): return (z[1][k]+1)*b<=a*2**z[0][k]
def cov(z): return (ge(z,0,7,8) and le(z,1,1,8)) or (ge(z,0,7,8) and ge(z,2,7,8)) or (le(z,1,1,8) and ge(z,2,7,8)) or (ge(z,0,3,4) and ge(z,1,3,4))
print(sum(1 for z in R if not cov(z)))

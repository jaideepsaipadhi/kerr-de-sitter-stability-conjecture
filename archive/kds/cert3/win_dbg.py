import pickle, numpy as np
exec(open('win_fac.py').read().split("for k,(F1,Y,P)")[0])
exec(open('/tmp/lz0c.py').read().split("for rho in")[0])
F1,Y,P=res['H2']; fP=fe(P); fF=fe(F1); fY=fe(Y)
pv=fP(X,1/V,T); idx=np.where(pv<0)[0][:5]
for i in idx:
    g=geom(X[i],V[i]); Q2,Q1,Q0,z0,A2,C0,t=coeffs(g,T[i])
    s=np.sqrt(fe({k:v for k,v in {}.items()}) if False else 0)
    print(X[i],V[i],T[i],'H2 direct/C0',(Q2+Q1/2+Q0/4)/C0,'P',pv[i],'F1',fF(X[i:i+1],1/V[i:i+1],T[i:i+1]),'Y',fY(X[i:i+1],1/V[i:i+1],T[i:i+1]))
print('neg count',np.sum(pv<0),'of',len(pv))

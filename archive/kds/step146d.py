import numpy as np, pickle
bad=pickle.load(open('bad146.pkl','rb'))
for P in bad[:8]:
    tt,zz,AA,BB,Sh=P; zh=zz-Sh; mu_=zh*(zh-1)*(tt-zh); mup=(zh-1)*(tt-zh)+zh*(tt-zh)-zh*(zh-1); H=2*mu_+Sh*mup; w0=H/Sh
    Lh=AA/zh+BB-3*zh
    u=np.linspace(0,1,2001)[:-1]; q=zh+Sh*u; s=zz-q; L=AA/q+BB-3*q; Lp=AA/q**2+3; mu=q*(q-1)*(tt-q)
    # pointwise admissible theta interval: g(th)=th(2-th)s2L2+th mu L-(2-th)mu s Lp >0
    a=-s*s*L*L; b=2*s*s*L*L+mu*L+mu*s*Lp; c=-2*mu*s*Lp   # a th^2 + b th + c
    disc=b*b-4*a*c; thlo=(-b+np.sqrt(disc))/(2*a); thhi=(-b-np.sqrt(disc))/(2*a)
    print('w0/(sL)hat %.3f'%(w0/(Sh*Lh)),'min thlo-needed max over z %.3f at u=%.3f'%(thlo.max(),u[thlo.argmax()]),'thhi min %.3f'%thhi.min(), 'u where thlo>w0/(sL):',(thlo>w0/(Sh*Lh)).any())

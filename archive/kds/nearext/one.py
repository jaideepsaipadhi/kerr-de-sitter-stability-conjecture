exec(open('/tmp/sweep.py').read().split("rng=")[0])
d=data(0.4545,1.0,5.66,1,0.986); t,A,B,z0,q1=d['t'],d['A'],d['B'],d['z0'],abs(d['q1'])
mu=lambda z: z*(z-1)*(t-z); mup=lambda z: (z-1)*(t-z)+z*(t-z)-z*(z-1)
L=lambda z: A/z+B-3*z; Lp=lambda z:-A/z**2-3; q=lambda z:q1*(z-z0)
def rhs(z,y):
    f,f1=y; return [f1,-((0.5*mup(z)*q1+q(z)*L(z))*f1+q(z)*Lp(z)*f)/(0.5*mu(z)*q1)]
p=2*(z0-1)*L(1)/(t-1); e=1e-6*(z0-1)
s=solve_ivp(rhs,[1+e,z0-1e-9],[e**p,p*e**(p-1)],rtol=1e-10,atol=1e-300)
print('p',p,'fmin>0',s.y[0].min()>0,'seam f1>=0',s.y[1][-1]>=0, 'x0',z0-1,'B',B,'floor',d['floor'])

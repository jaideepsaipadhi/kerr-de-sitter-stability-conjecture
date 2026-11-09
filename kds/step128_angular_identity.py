# Exact identity: angular KdS Teukolsky operator H (Suzuki-Takasugi-Umetsu eq. 2.3) = H0 + K_alpha + E, where
# H0 = -d(1-x^2)d + (m+sx)^2/(1-x^2) (spectrum l(l+1)-s^2), K_alpha = -d(al x^2(1-x^2))d >= 0, and
# E = N/(1+al x^2) - s al + 2 al x^2 with N as below (potential terms only).
import sympy as sp
x,al,xi,m,s=sp.symbols('x alpha xi m s')
VH=-(1+al)**2*xi**2/al+2*al*x**2-(1+al)/(1+al*x**2)*(2*s*(al*m-(1+al)*xi)*x-(1+al)**2*xi**2/al+2*m*(1+al)*xi+s**2) \
   +(1+al)**2*m**2/((1+al*x**2)*(1-x**2))+(1+al)*(s**2+2*s*m*x)/(1-x**2)
V0=(m+s*x)**2/(1-x**2)
N=-s**2*al*(1-x**2)+(1+al)**2*(xi**2*(1-x**2)+2*xi*(m+s*x))-2*s*al*m*(1+al)*x
E=N/(1+al*x**2)-s*al+2*al*x**2
print('VH - V0 - E =', sp.simplify(VH-V0-E))
sing=al*(m+s)**2/(2*(1-x))+al*(m-s)**2/(2*(1+x))
print('VH - V0 - sing - E =', sp.simplify(VH-V0-sing-E))

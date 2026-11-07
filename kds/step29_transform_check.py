"""
Step 29: check the Euler transform as a Riemann-Liouville integral based at z = 1.

Take the Frobenius solution of the original Heun equation at z = 1 with exponent rho = 1 - delta (the event-horizon QNM branch):
    y(w) = (w-1)^rho sum_k c_k (w-1)^k.
Transform termwise with mu = A - 2 (A the real-dual branch):
    ytilde(z) = sum_k c_k B(mu+1, rho+k+1) (z-1)^{mu+rho+1+k},
which is the Pochhammer-regularised RL integral int_1^z (z-w)^mu y(w) dw.
Checks:
  (1) ytilde solves the dual Heun equation (residual of the ODE at sample z);
  (2) its leading exponent is mu + rho + 1 = A - delta = 1 - delta_tilde, the dual's non-analytic branch;
  (3) Im(mu) != 0 and Im(rho) != 0 (no integer degeneracy), and the leading coefficient is nonzero.
"""
import warnings; warnings.filterwarnings('ignore')
import mpmath as mp
import sympy as sp
from step3_dual import teukolsky_heun, heun_params_from_operator, real_branch
from step5_window import ctdc_window

mp.mp.dps = 40


def frobenius_at_1(g, d, e, AB, q, t, K=80):
    """Series of the Heun equation about z = 1 with exponent rho = 1 - d.  Shift x = z - 1, so the singular points
    are at x = -1, 0, t - 1, and multiply the equation through by x(x+1)(x-(t-1))."""
    rho = 1 - d
    T = t - 1
    # Heun: (z)(z-1)(z-t) y'' + [g(z-1)(z-t) + d z(z-t) + e z(z-1)] y' + (AB z - q) y = 0, with z = 1 + x
    x = sp.symbols('x')
    P2 = sp.expand((1 + x)*x*(x - T))
    P1 = sp.expand(g*x*(x - T) + d*(1 + x)*(x - T) + e*(1 + x)*x)
    P0 = sp.expand(AB*(1 + x) - q)
    p2 = [complex(P2.coeff(x, j)) for j in range(4)]
    p1 = [complex(P1.coeff(x, j)) for j in range(3)]
    p0 = [complex(P0.coeff(x, j)) for j in range(2)]
    # y = sum c_k x^{k+rho}:  sum over j of coefficient * (n-j+rho)(...) c_{n-j}, collected at power x^{n+rho-1}
    c = [mp.mpc(1)]
    def F(j, nu):  # contribution of x^j terms acting on x^nu, mapped to power x^{nu+j-1}
        val = mp.mpc(0)
        if 0 <= j < 4: val += p2[j]*nu*(nu - 1) if j >= 0 else 0
        if 0 <= j - 1 + 1 < 3: val += p1[j]*nu if j < 3 else 0
        if 1 <= j + 1 - 1 and j - 1 < 2: val += p0[j - 1] if 0 <= j - 1 < 2 else 0
        return val
    # power balance: term x^j * (x^nu)'' -> x^{nu+j-2}; x^j (x^nu)' -> x^{nu+j-1}; x^j x^nu -> x^{nu+j}
    def coef(n_total):
        pass
    # direct linear recursion: collect power x^{rho+n-1}
    for n in range(1, K):
        acc = mp.mpc(0)
        for kk in range(0, n):
            nu = kk + rho
            sh = n - kk  # target power rho+n-1 = nu + sh - 1
            # x^j y'': power nu + j - 2 = nu + sh - 1 -> j = sh + 1
            j = sh + 1
            if 0 <= j < 4: acc += p2[j]*nu*(nu - 1)*c[kk]
            # x^j y': power nu + j - 1 -> j = sh
            j = sh
            if 0 <= j < 3: acc += p1[j]*nu*c[kk]
            # x^j y: power nu + j -> j = sh - 1
            j = sh - 1
            if 0 <= j < 2: acc += p0[j]*c[kk]
        nu = n + rho
        lead = p2[1]*nu*(nu - 1) + p1[0]*nu  # j = 1 for y'', j = 0 for y' (sh = 0)
        c.append(-acc/lead)
    # check the indicial equation for rho
    ind = p2[1]*rho*(rho - 1) + p1[0]*rho
    return rho, c, abs(ind)


def heun_residual(par, t, yfun, z, h=mp.mpf('1e-12')):
    g, d, e, AB, q = par
    y0 = yfun(z); y1 = (yfun(z + h) - yfun(z - h))/(2*h); y2 = (yfun(z + h) - 2*y0 + yfun(z - h))/h**2
    lhs = z*(z - 1)*(z - t)*y2 + (g*(z - 1)*(z - t) + d*z*(z - t) + e*z*(z - 1))*y1 + (AB*z - q)*y0
    scale = abs(z*(z - 1)*(z - t)*y2) + abs((g*(z - 1)*(z - t) + d*z*(z - t) + e*z*(z - 1))*y1) + abs((AB*z - q)*y0)
    return abs(lhs)/scale


if __name__ == '__main__':
    for (Lam, a, s, wf, m, lam) in [(0.05, 0.6, -2, 0.5, 2, 7.3), (0.02, 0.9, -2, 0.3, 3, 12.0), (0.05, 0.6, +2, 0.5, 2, 7.3)]:
        lo, hi, _, _ = ctdc_window(Lam, a, 1.0)
        om = m*(lo + wf*(hi - lo))
        data = teukolsky_heun(Lam, a, 1.0, s, complex(om), m, lam)
        e_ = data['exps']; t = mp.mpf(str(float(data['t'])))
        h = heun_params_from_operator(data['Pz'], data['Qz'], data['z'], data['t'], e_['0'][0], e_['1'][0], e_['t'][0],
                                      {'e': 2*s + 1, 'zinf': sp.Float(str(data['zinf']), 35)}, s)
        which, ht = real_branch(h, complex(t))
        A = ht['A']; mu = A - 2
        g, d, e, AB, q = [mp.mpc(h[k]) for k in ('g', 'd', 'e', 'AB', 'q')]
        rho, c, ind = frobenius_at_1(g, d, e, AB, q, t)
        radius = min(1, float(t) - 1)
        y = lambda zz: sum(c[k]*(zz - 1)**(k + rho) for k in range(len(c)))
        coeffs = [c[k]*mp.beta(mu + 1, rho + k + 1) for k in range(len(c))]
        yt = lambda zz: sum(coeffs[k]*(zz - 1)**(mu + rho + 1 + k) for k in range(len(c)))
        tpar = [mp.mpc(ht[k]) for k in ('g', 'd', 'e', 'AB', 'q')]
        zs = [1 + radius*f for f in (0.15, 0.3, 0.45)]
        res_orig = max(heun_residual((g, d, e, AB, q), t, y, mp.mpf(zz)) for zz in zs)
        res_dual = max(heun_residual(tpar, t, yt, mp.mpf(zz)) for zz in zs)
        lead_exp = mu + rho + 1; dual_branch = 1 - mp.mpc(ht['d'])
        print(f'Lam={Lam} a={a} s={s}: branch={"AB"[which]}  indicial residual={float(ind):.1e}  '
              f'ODE residual: original {float(res_orig):.1e}, dual {float(res_dual):.1e}')
        print(f'    leading exponent mu+rho+1 = {complex(lead_exp):.6f}  vs  dual 1-delta~ = {complex(dual_branch):.6f}')
        print(f'    Im mu = {float(mp.im(mu)):+.4f}, Im rho = {float(mp.im(rho)):+.4f}, |B(mu+1, rho+1)| = {float(abs(coeffs[0])):.3e}')

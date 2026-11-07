# Kerr–de Sitter spin-2 mode stability: dual + Hintz two-current pipeline

Requirements: `pip install sympy mpmath numpy scipy`

| script | what it does |
|---|---|
| step1_radial.py | Radial Teukolsky on KdS. Fixes the r^2-term convention (regularity at r = inf) and checks the Kerr limit |
| step2_euler.py | Derives every Euler-transform closure of the Heun operator (kernel (z-w)^{A-2} or (z-w)^{B-2}) |
| step3_dual.py | Teukolsky -> Heun (validated) -> all 16 gauge/transform duals; finds the ones with a REAL potential |
| step4_hintz_form.py | Writes the real duals as mu*W = qt^2 - Lt*mu |
| step5_window.py, step5b_bisect.py | Dual superradiance (qt(1) qt(t) < 0) vs the CTdC window (3.37): they agree to 10 digits |
| step6_hintz_lp.py | Hintz's two-current conditions on the dual, as LP feasibility over multipliers f |
| step7_verify.py | Re-checks the LP multipliers on a fine grid; `quick` or `full` sweep |
| step7b_diag.py, step7c_nearkerr.py | Failure diagnostics |
| step8_piecewise.py | Piecewise multipliers (continuity of f, kink condition q'[f'] <= 0), fine-grid check |
| step9_freeq.py | **Current method**: qt free in the interior (only the horizon values are fixed), L = (qt^2 - P4)/mu; optimise qt; cutting-plane fine verification |
| step9b_spot.py | Spot checks on extra geometries |

Run the current sweep:  `python3 step9_freeq.py medium > sweep_medium.txt`  (parallel; ~20-40 min)

## Explicit-structure track (steps 10-12)
| script | what it does |
|---|---|
| step10_structure.py | The straight-line qt already closes 10/11 hard cases; one extra constant closes 11/11 |
| step11_explicit.py | Cosmological side in closed form: qt linear, f+ = 1. Holds on all 672 grid cases |
| step12_event.py | Event side: f = u^p(1+p-pu), u=(z-1)/(z0-1). Closes 548/672 on the full grid and 760/840 with physical lambda = l(l+1)-4 |

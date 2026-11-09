#!/bin/sh
# Item 119: exact 2-D Bernstein certificate for (P1).  Charts cover [0,1]^2 in (T = 1/t, zeta):
python3 p1chart.py R1        # T in [0.1,1], zeta in [0,1]
python3 p1chart.py R2        # T in [0,0.1], zeta in [0.1,1]
python3 p1chart.py A 1/10    # corner, zeta <= T  (zeta = T w)
python3 p1chart.py B 1/10    # corner, T <= zeta  (T = zeta w)

#!/bin/sh
# Item 114: complete exact certificate of B3 on the A = 0 surface (bulk part), three charts, total < 10 s.
python3 chart2.py 0.05 0 0 0 1 1 1 1        # t <= 20, all zeta, y, B          (chart: bt = floor/B)
python3 chart1.py 0 0 0 0 0.05 0.9 1 1      # t >= 20, zeta <= 0.9              (chart: omega = 10/w)
python3 chart2.py 0 0.9 0 0 0.05 1 1 1      # t >= 20, zeta >= 0.9              (chart: bt)

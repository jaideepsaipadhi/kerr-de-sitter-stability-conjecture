#!/bin/sh
# LEMMAS 70: exact mult_(0,0) of the constraint minor delta_C >= 15 (K=4 Toeplitz on every line), three test points.
# Expected output: three lines "C: T_4 deficiency on every line = 15".  Roughly 1 h per point on one core; runs 2 in parallel.
cd "$(dirname "$0")"
python3 step59b_origin_mult.py 0 4 C > origin_C_pt0.log 2>&1 &
python3 step59b_origin_mult.py 1 4 C > origin_C_pt1.log 2>&1
wait
python3 step59b_origin_mult.py 2 4 C > origin_C_pt2.log 2>&1
grep -h "deficiency" origin_C_pt*.log

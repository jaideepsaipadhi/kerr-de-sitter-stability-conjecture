#!/bin/bash
# rerun 2D corner certificates with per-entry tolerance
export CONDS=q0,q1,v1,v2,c2neg,disc,mm,pm,mp
for prog in cover2d_f cover2d_f6; do for ch in EA EB; do
  for k in $(seq 1 8); do out=$(timeout 300 python3 $prog.py $ch 100 60 2>&1 | grep -i "closed\|open" | tail -1); echo "$prog $ch run$k: $out"
    echo "$out" | grep -q "open 0 " && break; done; done; done

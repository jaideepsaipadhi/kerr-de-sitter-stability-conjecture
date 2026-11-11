#!/bin/bash
while pgrep -f rerun_c2.sh >/dev/null; do sleep 5; done
for ch in "" KA KB; do for k in $(seq 1 6); do
  out=$(timeout 300 python3 cover3d_f.py margin/scn3_res.pkl Q2,H2,H1 scn3 100 60 $ch 2>&1 | grep -i "closed" | tail -1); echo "scn3 [$ch] run$k: $out"
  echo "$out" | grep -q "open 0 " && break; done; done

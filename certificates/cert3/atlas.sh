#!/bin/bash
# usage: atlas.sh res targets alt name
RES=$1; TG=$2; ALTS=$3; NM=$4
for c in ${CHARTS:-KA KB VT1 VT2 T1a T1b T2 T3a T3b XT1 XT2 U1 U2 U3a U3b EA EB}; do
  rm -f margin/c3f_${NM}${c}.json
  for k in 1 2 3 4 5 6; do
    if [ -z "$c" ]; then out=$(PERM2=1 NOEX=1 ALT=$ALTS timeout 400 python3 cover3d_h.py $RES $TG $NM 100 60 | grep closed);
    else out=$(PERM2=1 ALT=$ALTS timeout 400 python3 cover3d_h.py $RES $TG $NM 100 60 $c | grep closed); fi
    echo "$NM [$c] run$k: $out"
    if echo "$out" | grep -q "open 0 "; then break; fi
  done
done

#!/bin/bash
# usage: render_views.sh <prim> <outdir> [views...]  -- the make_solid_views view set, PNG + SVG
LV=$HOME/Projects/solid/solid_dd4hep/bin/solid_lineview
prim=$1; out=$2; shift 2
views=${@:-section_yz section_xz cutaway lines_yz}
mkdir -p $out
for v in $views; do
  case $v in
    section_yz) opts="--mode section --section x=0 --view side --fill --scale-bar";;
    section_xz) opts="--mode section --section y=0 --view top --fill --scale-bar";;
    cutaway)    opts="--mode hlr --dawncut -1 0 0 0 --view iso --perspective 40000 --no-traj --no-hits";;
    lines_yz)   opts="--mode hlr --dawncut -1 0 0 0 --view side";;
    beam_*)     z=${v#beam_}; opts="--mode section --section z=${z}cm --view beam --fill --scale-bar";;
  esac
  s=$(date +%s.%N)
  nice -n 10 python3 $LV $prim $opts -o $out/$v.png -o $out/$v.svg -q > $out/$v.log 2>&1
  echo "$v rc=$? $(echo "$(date +%s.%N) - $s" | bc) s"
done

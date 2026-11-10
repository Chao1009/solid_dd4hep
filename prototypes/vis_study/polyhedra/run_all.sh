#!/bin/bash
# G4 (standalone) and TGeo (fixed) tessellation of all 7 configurations, sequential, nice'd.
W=$SCRATCH/vis_study/polyhedra
I=$SCRATCH/verify_clone/install
P=$SCRATCH/overview_job/prims
SIF=$HOME/Projects/solid/containers/eic_xl-26.10.0-stable.sif
cd $W
for c in solid_SIDIS_He3 solid_SIDIS_NH3 solid_PVDIS solid_JPsi solid_DDVCS solid_JPsi_lite solid_DDVCS_lite; do
  for tool in g4poly tgeopoly; do
    pre=g4; [ $tool = tgeopoly ] && pre=tgeo
    s=$(date +%s.%N)
    nice -n 10 singularity exec $SIF bash -c ". /etc/profile >/dev/null 2>&1; source $I/bin/thissolid.sh $c >/dev/null; cd $W; python3 $tool.py -c \$DETECTOR_PATH/$c.xml -o out/all/${pre}_$c --bbox-from $P/$c/g4_0000.prim" > out/all/${pre}_$c.log 2>&1
    echo "$c $tool rc=$? wall=$(echo "$(date +%s.%N) - $s" | bc)"
  done
  nice -n 10 python3 compare_prims.py $P/$c/g4_0000.prim out/all/g4_$c.prim > out/all/cmp_dawn_g4_$c.txt 2>&1
  nice -n 10 python3 compare_npz.py out/all/g4_$c.npz out/all/tgeo_$c.npz > out/all/cmp_g4_tgeo_$c.txt 2>&1
  nice -n 10 python3 topo_stats.py out/all/g4_$c.npz out/all/tgeo_$c.npz > out/all/topo_$c.txt 2>&1
done
echo ALLDONE

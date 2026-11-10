#!/bin/bash
# usage: run_batch.sh RUNDIR OUTFILE  (env: STEER, CONFIG, SEED, GUN_*, NEVT)
SIF=$HOME/Projects/solid/containers/eic_xl-26.10.0-stable.sif
INST=$SCRATCH/verify_clone/install
CONFIG=${CONFIG:-solid_JPsi_lite}
RUNDIR=$(readlink -f "$1"); OUT=$2
mkdir -p "$RUNDIR"; cd "$RUNDIR" || exit 1
singularity exec "$SIF" bash -c ". /etc/profile >/dev/null 2>&1; unset DISPLAY WAYLAND_DISPLAY XAUTHORITY; \
  source $INST/bin/thissolid.sh $CONFIG; export SOLID_OPTICAL=0; cd $RUNDIR; \
  nice -n 10 ddsim --steeringFile ${STEER} --compactFile \$DETECTOR_PATH/\$DETECTOR_CONFIG.xml --runType batch \
     --outputFile $RUNDIR/$OUT > $RUNDIR/ddsim_batch.log 2>&1; echo rc=\$?"

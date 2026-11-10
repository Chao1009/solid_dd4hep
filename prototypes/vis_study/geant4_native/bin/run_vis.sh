#!/bin/bash
# usage: run_vis.sh RUNDIR MACRO [extra ddsim args...]
# Runs ddsim --runType vis headless (DISPLAY unset) in the eic_xl container with SoLID-Lite J/psi.
SIF=$HOME/Projects/solid/containers/eic_xl-26.10.0-stable.sif
INST=$SCRATCH/verify_clone/install
CONFIG=${CONFIG:-solid_JPsi_lite}
RUNDIR=$(readlink -f "$1"); MAC=$(readlink -f "$2"); shift 2
mkdir -p "$RUNDIR"
cd "$RUNDIR" || exit 1
EXTRA="$*"
singularity exec "$SIF" bash -c ". /etc/profile >/dev/null 2>&1; unset DISPLAY WAYLAND_DISPLAY XAUTHORITY; \
  source $INST/bin/thissolid.sh $CONFIG; export SOLID_OPTICAL=0 G4DAWNFILE_VIEWER=true; cd $RUNDIR; \
  start=\$(date +%s.%N); \
  echo exit | nice -n 10 ddsim --steeringFile ${STEER:-\$DETECTOR_PATH/steering/solid_sim.py} \
     --compactFile \$DETECTOR_PATH/\$DETECTOR_CONFIG.xml --runType vis --macroFile $MAC \
     --outputFile $RUNDIR/sim.edm4hep.root $EXTRA > $RUNDIR/ddsim.log 2>&1; rc=\$?; \
  end=\$(date +%s.%N); echo \"rc=\$rc wall=\$(python3 -c \"print(round(\$end-\$start,1))\") s\" | tee -a $RUNDIR/ddsim.log"

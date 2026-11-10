#!/bin/bash
# usage: run_reve.sh EVENT_TXT LOG SECONDS [extra root flags]   (runs inside the container)
source $SCRATCH/vis_study/web_interactive/scripts/env.sh
unset DISPLAY
cd $W/reve2
root -l $4 -q "$W/scripts/solid_reve.C(\"$W/geo2/solid_JPsi_lite_display.root\",\"$1\",8795,$3)" > $2 2>&1

# common environment for the web_interactive prototypes (source inside the container)
. /etc/profile >/dev/null 2>&1
export W=$SCRATCH/vis_study/web_interactive
export INST=$SCRATCH/verify_clone/install
source $INST/bin/thissolid.sh solid_JPsi_lite >/dev/null
export SOLID_OPTICAL=0

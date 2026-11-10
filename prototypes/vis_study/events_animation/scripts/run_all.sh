#!/bin/bash
# Reproduce the events_animation prototype end to end, with logs and timings in logs/.
#   run_all.sh [stage ...]    stages: plugin sim bench extract render animate encode svgcheck (default: all)
# Container stages run in eic_xl 26.10.0-stable; rendering/encoding on the host (python3 + numpy + pycairo +
# Pillow, GStreamer 1.20, Firefox for the SVG check).
W=$SCRATCH/vis_study/events_animation
SIF=$HOME/Projects/solid/containers/eic_xl-26.10.0-stable.sif
INST=$SCRATCH/verify_clone/install
PRIM=$SCRATCH/overview_job/prims/solid_JPsi_lite/g4_0000.prim
STAGES=${*:-plugin sim bench extract render animate encode svgcheck}
L=$W/logs; mkdir -p $L $W/out
TIME="/usr/bin/time -f %e_s_wall_%M_kB_maxRSS"
C() { singularity exec $SIF bash -c ". /etc/profile >/dev/null 2>&1; source $INST/bin/thissolid.sh solid_JPsi_lite >/dev/null; export LD_LIBRARY_PATH=$W/plugin:\$LD_LIBRARY_PATH SOLID_OPTICAL=0; $1"; }
REC='{"name":"SolidStepRecorder/StepRec","parameter":{"OutputFile":"steps.bin","MaxDistance":20.0,"MaxAngle":0.04}}'
DDSIM="nice -n 10 $TIME ddsim --steeringFile \$DETECTOR_PATH/steering/solid_sim.py --compactFile \$DETECTOR_PATH/\$DETECTOR_CONFIG.xml"
EV="env PYTHONDONTWRITEBYTECODE=1 nice -n 10 $TIME python3 $W/scripts/evdisp.py --prim $PRIM --geo-cache $W/sim/JPsi_lite_geo_v03.pkl"
for st in $STAGES; do case $st in
plugin)   # 1. DDG4 stepping action -> libSolidStepRecorder.so + .components
  C "cd $W/plugin && $TIME g++ -O2 -fPIC -shared SolidStepRecorder.cpp -o libSolidStepRecorder.so \$(root-config --cflags) \$(geant4-config --cflags) -I/opt/local/include -L/opt/local/lib -lDDG4 -lDDCore \$(geant4-config --libs) \$(root-config --libs) && listcomponents_dd4hep -o libSolidStepRecorder.components libSolidStepRecorder.so && cat libSolidStepRecorder.components" > $L/plugin.log 2>&1 ;;
sim)      # 2. ddsim + recorder: 5 x 4 GeV e- (default gun of solid_sim.py), 2 hand-made J/psi-like events (HepMC3)
  C "mkdir -p $W/sim/egun5 && cd $W/sim/egun5 && $DDSIM -N 5 --random.seed 12345 --outputFile egun5.edm4hep.root --action.step '$REC'" > $L/sim_egun5.log 2>&1
  C "mkdir -p $W/sim/jpsi_like && cd $W/sim/jpsi_like && cp $W/inputs/jpsi_like.hepmc . && SOLID_GUN=0 $DDSIM -N 2 --random.seed 777 --inputFiles jpsi_like.hepmc --outputFile jpsi_like.edm4hep.root --action.step '$REC'" > $L/sim_jpsi_like.log 2>&1 ;;
bench)    # 2b. recorder overhead: same 5 events without / with a coarser recorder
  C "mkdir -p $W/sim/bench && cd $W/sim/bench && $DDSIM -N 5 --random.seed 12345 --outputFile bench_none.edm4hep.root" > $L/bench_none.log 2>&1
  C "mkdir -p $W/sim/bench && cd $W/sim/bench && $DDSIM -N 5 --random.seed 12345 --outputFile bench_coarse.edm4hep.root --action.step '{\"name\":\"SolidStepRecorder/StepRec\",\"parameter\":{\"OutputFile\":\"steps_coarse.bin\",\"MaxDistance\":50.0,\"MaxAngle\":0.1,\"MinEkin\":1.0}}'" > $L/bench_coarse.log 2>&1 ;;
extract)  # 3. EDM4hep + DD4hep (IDDescriptor, CellIDPositionConverter, VolumeManager, placements) -> npz
  C "cd $W/sim/egun5 && $TIME python3 $W/scripts/extract_event.py \$DETECTOR_PATH/\$DETECTOR_CONFIG.xml egun5.edm4hep.root egun5_event.npz" > $L/extract_egun5.log 2>&1
  C "cd $W/sim/jpsi_like && $TIME python3 $W/scripts/extract_event.py \$DETECTOR_PATH/\$DETECTOR_CONFIG.xml jpsi_like.edm4hep.root jpsi_like_event.npz" > $L/extract_jpsi_like.log 2>&1 ;;
render)   # 4. static figures
  for e in 0 2; do $EV --event-npz $W/sim/egun5/egun5_event.npz --steps $W/sim/egun5/steps.bin -e $e -o $W/out --tag ev$e \
      --views side,sidezoom,beam,cutaway --formats svg,png,pdf > $L/render_ev$e.log 2>&1; done
  for e in 0 1; do $EV --event-npz $W/sim/jpsi_like/jpsi_like_event.npz --steps $W/sim/jpsi_like/steps.bin -e $e -o $W/out \
      --tag jpsi$e --views side,beam,cutaway --formats svg,png,pdf > $L/render_jpsi$e.log 2>&1; done
  $EV --event-npz $W/sim/jpsi_like/jpsi_like_event.npz --steps $W/sim/jpsi_like/steps.bin -e 1 -o $W/out --tag jpsi1_ortho \
      --proj ortho --views side --formats png > $L/render_jpsi1_ortho.log 2>&1
  $EV --event-npz $W/sim/egun5/egun5_event.npz --traj edm4hep -e 2 -o $W/out --tag ev2_routeB --views side,sidezoom \
      --formats png > $L/render_ev2_routeB.log 2>&1
  $EV --event-npz $W/sim/jpsi_like/jpsi_like_event.npz --traj edm4hep -e 1 -o $W/out --tag jpsi1_routeB --views side \
      --formats png > $L/render_jpsi1_routeB.log 2>&1 ;;
animate)  # 5. PNG frames (cached geometry layers) + animated SVG
  rm -rf $W/out/frames_*
  $EV --event-npz $W/sim/egun5/egun5_event.npz --steps $W/sim/egun5/steps.bin -e 2 -o $W/out --tag ev2 --views side,beam,cutaway \
      --formats png --animate side,beam,cutaway --frames 120 --frame-scale 0.6 --svg-anim > $L/animate_ev2.log 2>&1
  $EV --event-npz $W/sim/jpsi_like/jpsi_like_event.npz --steps $W/sim/jpsi_like/steps.bin -e 1 -o $W/out --tag jpsi1 \
      --views side,beam --formats png --animate side,beam --frames 120 --frame-scale 0.6 --svg-anim > $L/animate_jpsi1.log 2>&1 ;;
encode)   # 6. GIF / WebP (Pillow), WebM (GStreamer vp9enc); MP4 attempted (needs an H.264/AV1/MPEG-4 encoder)
  cd $W/out && for v in ev2_side ev2_cutaway ev2_beam jpsi1_side jpsi1_beam; do nice -n 10 python3 $W/scripts/encode.py frames_$v ${v}_anim \
      --fps 12 --formats gif,webp,webm,mp4 --webp-lossless; done > $L/encode.log 2>&1 ;;
svgcheck) # 7. animated SVG rendered by headless Firefox at fixed animation times (SMIL paused + setCurrentTime)
  nice -n 10 python3 $W/scripts/svgcheck.py $W/out/ev2_side_anim.svg $W/out/svgcheck 9.5 26.1 26.8 28.2 > $L/svgcheck.log 2>&1
  nice -n 10 python3 $W/scripts/svgcheck.py $W/out/jpsi1_side_anim.svg $W/out/svgcheck 5 12 22 29 > $L/svgcheck_jpsi1.log 2>&1 ;;
esac; done

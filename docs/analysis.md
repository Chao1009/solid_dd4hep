# Analysis of simulated events: particle tables and acceptance

The first layer between the simulation and the physics projections. It works on simulated hits with the Monte
Carlo truth, not on reconstructed quantities.

| step | tool | output |
|---|---|---|
| single particles, flat in momentum, polar angle, azimuth and along the target | `scripts/analysis/gun_events.py` | HepMC3 file |
| simulation | `ddsim` with `steering/solid_sim.py` | EDM4hep file |
| one row per generated particle with what it left in every detector | `scripts/analysis/particle_table.py` | `table.npz` |
| acceptance in (θ, p), with or without azimuthal segments, compared with GEMC | `scripts/analysis/acceptance.py` | maps, figure |

```
gun_events.py --pdg 11 -n 10000 --p 0.5:8 --theta 6:28 --z=-370:-330 -o e.hepmc
SOLID_GUN=0 SOLID_OPTICAL=0 ddsim --steeringFile $DETECTOR_PATH/steering/solid_sim.py --inputFiles e.hepmc -N 10000 \
        --outputFile sim.edm4hep.root
particle_table.py sim.edm4hep.root -o table.npz
acceptance.py table.npz --config full --plot acceptance.png --gemc acceptance_solid_SIDIS_He3_electron_1e7_201701_output_final.root
```

This is the method of the GEMC studies (`solid_gemc/analysis/acceptance`): particles flat in P and θ, even along
the target, and a map of the detected fraction. The maps of both simulations can therefore be compared bin by
bin, and physics events can be folded with them in the same way as before.

## The particle table

`python/solid_ana/table.py`. For every generator particle (generatorStatus 1):

- kinematics at the vertex: `pdg`, `p`, `theta`, `phi`, `vx`, `vy`, `vz` (GeV, rad, mm);
- for every tracker-type collection (GEM planes, scintillator pads, MRPC) and layer: whether the particle left a
  hit, the position and time of its first hit, and the azimuthal element of that hit (the GEM sector, the
  scintillator pad);
- for every calorimeter collection: the visible energy from the particle and all its descendants, the position
  of the cell with most of it, and the energy of the whole event.

Two properties of ddsim matter here. Showers that start inside the tracking region, which holds the large-angle
calorimeter, are stored as particles of their own: the table adds the energy of all descendants. And the hits of
secondaries that are not stored are booked on their stored ancestor: shower leakage can therefore add tracker
hits behind a calorimeter to a particle's row.

## What counts as detected

`python/solid_ana/accept.py`. A particle is detected in an arm when it left hits in enough GEM planes of the arm
and visible energy in the arm's calorimeter; a scintillator pad can be asked for as well.

| arm | GEM planes | default requirement |
|---|---|---|
| forward, SoLID-Lite | 5 forward planes | 4 planes, more than 50 MeV visible in the forward shower modules |
| large angle, SoLID-Lite | 4 large-angle planes | 3 planes, more than 50 MeV visible in the large-angle shower modules |
| forward, full SoLID | planes 2–6 of the six | the same |
| large angle, full SoLID | planes 1–4 | the same |

The shashlik modules show about 20 % of an electron's energy, so 50 MeV visible is about 250 MeV. These are loose
requirements: they give acceptance, not efficiency. No Cherenkov is asked for, since optical photons are off in
these runs.

## Azimuthal segments as a selection

For SoLID-Lite the forward detectors may cover only segments of the azimuth
([lite_options.md](lite_options.md)). `acceptance.py --fa-segments 4x24` keeps a particle when the GEM sectors
it hit, its scintillator pad and its calorimeter cell lie inside the segments, which is what a detector built
only there would record. One simulation of a full-azimuth detector then serves every layout of segments. The
table stores detector elements, not the particle's azimuth at the vertex, because the solenoid turns positive
and negative particles in opposite directions on their way to the detectors.

This selection is exact as long as the removed elements were not needed for the particle to reach the kept ones.
The segmented geometry (`solid_SIDIS_He3_lite.xml`) exists to check that.

## Trapped particles

A low-energy electron born in a vacuum volume inside the solenoid field can be trapped between two regions of
stronger field. It loses no energy, every step completes, and Geant4's looper protection does not act: one such
event stopped a production job for good. `steering/solid_sim.py` therefore adds the `SoLIDTrackTimeLimit`
stepping action, which stops charged particles still being tracked 100 µs after the event
(`SOLID_MAX_TIME=<ns>`, 0 switches it off). The log reports each stopped track and their number.

## Test

ctest `analysis` (`scripts/analysis/test_analysis.py`): 80 electrons through `ddsim`; the table against the
generated kinematics; GEM sectors and scintillator pads against the hit positions; the segment selection against
the azimuth of the hits; the acceptance tool; the time limit.

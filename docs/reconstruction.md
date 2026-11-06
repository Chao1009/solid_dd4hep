# SoLID reconstruction framework and software environment

Recommendation of November 2026. The full report, with the inspection of the
`eic_xl:26.10.0-stable` image and the sources, was produced alongside this repository update.

## Recommendation

**Run EICrecon (JANA2) as the production reconstruction framework, as an out-of-tree set of SoLID plugins. Write
every SoLID-specific algorithm against the framework-neutral `eic/algorithms` interface, on EDM4hep/EDM4eic
types.** A Gaudi (k4FWCore) wrapper then costs about 40 lines per algorithm and can be added at any time. Build the
Gaudi side only if a real consumer and an owner appear (decision gate G2, mid 2027).

In short: "use EICrecon, but code so Gaudi stays a cheap option". The reasons, checked on 2026-11-03:

1. **Gaudi is no longer in the EIC software stack.** eic/containers PR #377 ("rm juggler/gaudi", merged
   2026-08-04) removed Gaudi and Juggler as "not used anymore by anything", and `eic/algorithms` v2.0.0
   (2026-05-12) removed its Gaudi/Juggler implementations. The `eic_xl:26.10.0-stable` image has no Gaudi,
   k4FWCore or Juggler.
2. **The algorithms SoLID needs already ship in EICrecon 1.42/1.43**, and their libraries do not depend on JANA, so
   any framework can reuse them:
   - calorimeter digitization and island clustering (hexagonal adjacency included);
   - MPGD strip digitization;
   - photomultiplier digitization for the gas Cherenkovs;
   - ACTS seeding, CKF, ambiguity solving and propagation;
   - PID lookup and ONNX inference.
3. **SoLID has already taken this path.**
   - The 2022 white paper asks for software "closely aligned with that of the EIC".
   - The January 2025 collaboration plan set the toolchain to EDM4hep + DD4hep + ACTS + EICrecon.
   - `solid_recon` (2025-26) already runs FAEC/LAEC clustering, GEM digitization and ACTS CKF tracking in the
     EICrecon structure. It is a vendored copy of EICrecon 25.05, however, and must be re-based as out-of-tree
     plugins.
4. **About 1 FTE of staffing** is enough to configure and extend one ecosystem. It is not enough to rebuild a chain
   in Key4hep or to keep two first-class framework integrations alive.
5. **"Gaudi is more mature" holds for the event-loop kernel, not for the SoLID-relevant reconstruction content.**
   Key4hep tracking (k4ActsTracking) needs hand-written per-detector ACTS blueprint code that SoLID would own
   alone. JANA2 is JLab's framework and brings the EVIO/CCDB/streaming heritage (GlueX) that real SoLID data will
   need.

| | A: Gaudi/Key4hep native | B: EICrecon out-of-tree plugins | C: neutral algorithms + JANA2 and Gaudi wrappers |
|---|---|---|---|
| effort to a first full chain | ~9-12 person-months | ~3-4 person-months (re-base solid_recon) | B + 4-6 person-weeks for the Gaudi side |
| reuse | little SoLID-like reconstruction | highest: generic factories, ACTS path matching the solid_dd4hep annotations | same as B |
| risk | SoLID would be the only NP/JLab Key4hep user; Gaudi was dropped from eic_xl | ePIC API churn (pin stable images, run nightly CI) | two configuration surfaces; thread-safety differences |
| container | Key4hep CVMFS, or eic_xl + 12 spack packages | eic_xl as is | eic_xl + Gaudi layer |
| verdict | not recommended | **production path** | **coding discipline**; build the Gaudi side after gate G2 |

What would change this: SoLID aligning with Key4hep/FCC instead of the EIC, JANA2 support ending, or EICrecon
governance blocking changes SoLID needs.

## Architecture (solid_recon, to be re-based)

```
solid_recon/
  algorithms/   framework-neutral SoLID algorithms (algorithms::Algorithm, EDM4hep/EDM4eic), Catch2 tests
  jana/         JOmniFactory wrappers + SoLID plugins (solid_ecal, solid_gem, solid_tracking, solid_cherenkov, ...)
  gaudi/        optional k4FWCore functional wrappers (CMake option SOLID_WITH_GAUDI, off by default)
  configs/      per-configuration JANA parameter files (sidis_he3.cfg, pvdis.cfg, jpsi_lite.cfg, ...)
```

| task | reuse from EICrecon | new SoLID code | phase |
|---|---|---|---|
| EC digitization and clustering | CalorimeterHitDigi/Reco, CalorimeterIslandCluster (hexagonal adjacency), ClusterRecoCoG | PreshowerMatch, per-configuration constants | 1 |
| GEM | MPGDTrackerDigi (v0), TrackerHitReconstruction | GEM u/v strip digitization (libsolgem physics), u/v hit pairing | 1-2 |
| tracking | CKFTracking, AmbiguitySolver, TrackPropagation | SoLIDTrackSeeding (from solid_recon), CKF configs, beamline vertex | 1-2 |
| LGC/HGC | PhotoMultiplierHitDigi | sector photoelectron sums, track matching | 1 |
| MRPC, SPD | CalorimeterHitDigi, pulse/CFD tools | MRPC and SPD digitization | 2 |
| trigger, PID, kinematics | MergeParticleID, ONNXInference | trigger emulator, PID combiner, fixed-target kinematics | 2 |
| real data | JANA2 event sources | EVIO/CODA sources, CCDB | 3 |

Phases (about 1 FTE, starting Nov 2026):

| phase | window | content |
|---|---|---|
| 0 | Nov-Dec 2026 | environment, CI, per-configuration compacts (this repository); regenerate the ACTS material map with ACTS 46 |
| 1 | Jan-Mar 2027 | first full chain, benchmarked against preCDR Table 22 |
| gate G1 | end of Mar 2027 | performance reproduces the GEMC-era results |
| 2 | Apr-Sep 2027 | realism: full GEM digitization, background overlay, MRPC/SPD, trigger, PID |
| gate G2 | mid 2027 | Gaudi decision |
| 3 | 2028+ | data readiness |

## Software environment ("VM")

eic-shell is an Apptainer/Docker container plus a launcher, not a VM. **A thin SoLID layer on top of eic_xl is worth
it; a custom software stack is not.**

* `containers/solid-shell` launches the eic_xl image (local `.sif`, CVMFS or docker) with this repository's install
  sourced (`thissolid.sh`). It works today without building anything.
* `containers/Dockerfile` builds `solid_xl` `FROM eicweb/eic_xl:<tag>-stable` and adds solid_dd4hep (and later
  solid_recon). A `WITH_GAUDI` build argument adds Gaudi + k4FWCore from the image's own spack stack; 12 extra
  packages were verified to concretize.
* Publishing: push to GHCR or Docker Hub, then add one line to the OSG cvmfs-singularity-sync list. On ifarm the
  image then appears under `/cvmfs/singularity.opensciencegrid.org/...`, and
  `eic-shell -o jeffersonlab -c solid_xl -v <tag>` works there.
* Effort: 3-5 person-days to set up, then about 0.5 day per month to bump the EIC tag.

The `.gitlab-ci.yml` of this repository targeted a dead Juggler image on eicweb. It is replaced by a GitHub Actions
workflow (`.github/workflows/ci.yml`) running in `eicweb/eic_xl`.
